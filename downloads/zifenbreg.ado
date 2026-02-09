*! zifenbreg v1.0.0 — Zero-inflated fixed-effects negative binomial regression
*! Concentrated/profile likelihood NB FE with logit zero-inflation
*! No FE in inflate equation (IPP constraint); FE via concentration in count

program define zifenbreg, eclass sortpreserve
    version 16.0

    if replay() {
        if "`e(cmd)'" != "zifenbreg" {
            error 301
        }
        Display `0'
        exit
    }

    Estimate `0'
end

program define Estimate, eclass

    syntax varlist(numeric fv ts min=2) [if] [in] , ///
        INFlate(varlist numeric fv ts)    ///
        [ NB1                             ///
          vce(string)                     ///
          Level(cilevel)                  ///
          TOLerance(real 1e-8)            ///
          ITERate(integer 1000)           ///
          NOLOG                           ///
          IRR                             ///
          EXPosure(varname numeric)       ///
          OFFset(varname numeric)         ///
          VUONG                           ///
        ]

    * --- Panel setup ---
    _xt, trequired
    local ivar `r(ivar)'
    local tvar `r(tvar)'

    * --- Parse depvar and indepvars ---
    gettoken depvar indepvars : varlist
    _fv_check_depvar `depvar'

    * --- Mark sample ---
    marksample touse
    markout `touse' `depvar' `indepvars' `inflate' `exposure' `offset'
    markout `touse' `ivar' `tvar', strok

    * --- Check depvar is nonneg integer ---
    capture assert `depvar' == int(`depvar') & `depvar' >= 0 if `touse'
    if _rc {
        di as error "`depvar' must be a nonnegative integer"
        exit 459
    }

    * --- Expand factor variables, remove collinear ---
    _rmcoll `indepvars' if `touse', noconstant expand
    local indepvars `r(varlist)'
    local k_count : word count `indepvars'

    if `k_count' == 0 {
        di as error "no regressors specified"
        exit 102
    }

    * --- Inflate variables (with constant) ---
    _rmcoll `inflate' if `touse', expand
    local inflate `r(varlist)'
    local k_inflate : word count `inflate'

    * k_inflate counts covariates only; +1 for constant added in Mata
    local k_inflate_total = `k_inflate' + 1

    * --- Model type ---
    local nb1_flag = ("`nb1'" != "")
    local modeltype = cond(`nb1_flag', "NB1", "NB2")

    * --- VCE parsing ---
    local vcetype "cluster"
    local clustvar "`ivar'"
    if `"`vce'"' != "" {
        gettoken vce1 vce2 : vce, parse(" ,")
        if "`vce1'" == "robust" {
            local vcetype "cluster"
            local clustvar "`ivar'"
        }
        else if "`vce1'" == "cluster" {
            local vcetype "cluster"
            local clustvar `vce2'
            if "`clustvar'" == "" {
                local clustvar "`ivar'"
            }
            confirm variable `clustvar'
        }
        else {
            di as error "vce() must be robust or cluster clustvar"
            exit 198
        }
    }

    * --- Offset / exposure ---
    tempvar offvar
    if "`exposure'" != "" & "`offset'" != "" {
        di as error "only one of exposure() or offset() allowed"
        exit 198
    }
    if "`exposure'" != "" {
        qui gen double `offvar' = ln(`exposure') if `touse'
        local offopt "offset(`exposure')"
        local offlab "ln(`exposure')"
    }
    else if "`offset'" != "" {
        qui gen double `offvar' = `offset' if `touse'
        local offopt "offset(`offset')"
        local offlab "`offset'"
    }
    else {
        qui gen double `offvar' = 0 if `touse'
    }

    * --- Drop singletons and all-zero groups ---
    tempvar Ti ysum
    sort `ivar'
    qui by `ivar': egen `Ti' = count(`depvar') if `touse'
    qui by `ivar': egen `ysum' = total(`depvar') if `touse'

    qui count if `touse' & `Ti' <= 1
    local n_single = r(N)
    if `n_single' > 0 {
        di as text "note: `n_single' singleton observations dropped"
        qui replace `touse' = 0 if `Ti' <= 1
    }

    qui count if `touse' & `ysum' == 0
    local n_allzero = r(N)
    if `n_allzero' > 0 {
        di as text "note: `n_allzero' observations in all-zero groups dropped"
        qui replace `touse' = 0 if `ysum' == 0
    }

    * --- Final sample counts ---
    qui count if `touse'
    local N = r(N)
    if `N' == 0 {
        di as error "no observations"
        exit 2000
    }

    tempvar grp
    qui egen `grp' = group(`ivar') if `touse'
    qui summ `grp' if `touse', meanonly
    local N_g = r(max)

    * --- Verbose? ---
    local verbose = ("`nolog'" == "")

    * --- Starting values for beta ---
    tempname b0_count
    local lnalpha_init 0
    capture qui poisson `depvar' `indepvars' if `touse', ///
        noconstant offset(`offvar') iterate(50)
    if _rc == 0 {
        matrix `b0_count' = e(b)
    }
    else {
        matrix `b0_count' = J(1, `k_count', 0)
    }

    * For NB1: get lnalpha start from pooled nbreg
    if `nb1_flag' {
        capture qui nbreg `depvar' `indepvars' if `touse', ///
            dispersion(constant) iterate(50)
        if _rc == 0 {
            capture local lnalpha_init = _b[/lnalpha]
            if missing(`lnalpha_init') local lnalpha_init 0
        }
    }

    * --- Starting values: logit for gamma ---
    tempvar yzero
    tempname b0_inflate
    qui gen byte `yzero' = (`depvar' == 0) if `touse'
    capture qui logit `yzero' `inflate' if `touse', iterate(50)
    if _rc == 0 {
        matrix `b0_inflate' = e(b)
    }
    else {
        matrix `b0_inflate' = J(1, `k_inflate_total', 0)
    }

    * --- Call Mata estimator ---
    tempname b_count b_inflate V lnalpha ll ll_0 converged ll_fenbreg

    mata: _zifenbreg_estimate(               ///
        "`depvar'",                           ///
        "`indepvars'",                        ///
        "`inflate'",                          ///
        "`offvar'",                           ///
        "`ivar'",                             ///
        "`grp'",                              ///
        "`touse'",                            ///
        `nb1_flag',                           ///
        "`clustvar'",                         ///
        `tolerance',                          ///
        `iterate',                            ///
        `verbose',                            ///
        "`b0_count'",                         ///
        "`b0_inflate'",                       ///
        `lnalpha_init',                       ///
        "`b_count'", "`b_inflate'", "`V'",    ///
        "`lnalpha'", "`ll'", "`ll_0'",       ///
        "`converged'", "`ll_fenbreg'"         ///
    )

    if `converged' == 0 {
        di as error "convergence not achieved"
        exit 430
    }

    * --- Build combined coefficient vector and VCE ---
    * Layout: [count:beta | inflate:gamma | lnalpha:_cons]
    local k_total = `k_count' + `k_inflate_total' + 1

    tempname b_full V_full
    matrix `b_full' = `b_count', `b_inflate', `lnalpha'

    * Column names: count eq, inflate eq, lnalpha eq
    local cnames_count ""
    foreach v of local indepvars {
        local cnames_count "`cnames_count' count:`v'"
    }
    local cnames_inflate ""
    foreach v of local inflate {
        local cnames_inflate "`cnames_inflate' inflate:`v'"
    }
    local cnames_inflate "`cnames_inflate' inflate:_cons"
    local cnames_lna "lnalpha:_cons"

    matrix colnames `b_full' = `cnames_count' `cnames_inflate' `cnames_lna'
    matrix rownames `b_full' = `depvar'

    matrix `V_full' = `V'
    matrix colnames `V_full' = `cnames_count' `cnames_inflate' `cnames_lna'
    matrix rownames `V_full' = `cnames_count' `cnames_inflate' `cnames_lna'

    * --- Post results ---
    ereturn post `b_full' `V_full', esample(`touse') depname(`depvar') obs(`N')

    ereturn local cmd "zifenbreg"
    ereturn local cmdline "zifenbreg `0'"
    ereturn local model "`modeltype'"
    ereturn local depvar "`depvar'"
    ereturn local ivar "`ivar'"
    ereturn local tvar "`tvar'"
    ereturn local vce "`vcetype'"
    ereturn local vcetype "Sandwich"
    ereturn local clustvar "`clustvar'"
    ereturn local predict "zifenbreg_p"
    ereturn local inflate "`inflate'"
    if "`exposure'" != "" ereturn local exposure "`exposure'"
    if "`offset'" != ""   ereturn local offset "`offset'"
    if "`offlab'" != ""   ereturn local offset1 "`offlab'"

    ereturn scalar N = `N'
    ereturn scalar N_g = `N_g'
    ereturn scalar k_count = `k_count'
    ereturn scalar k_inflate = `k_inflate_total'
    ereturn scalar ll = `ll'
    ereturn scalar ll_0 = `ll_0'
    ereturn scalar lnalpha = `lnalpha'
    ereturn scalar alpha = exp(`lnalpha')
    ereturn scalar converged = `converged'
    ereturn scalar k_eq = 3
    ereturn scalar k_dv = 1
    ereturn scalar rank = `k_total'
    ereturn scalar N_clust = `N_g'

    ereturn local title "Zero-inflated FE negative binomial regression (`modeltype')"
    ereturn local chi2type "Wald"
    tempname bw Vw chi2val
    matrix `bw' = e(b)
    matrix `bw' = `bw'[1, 1..`k_count']
    matrix `Vw' = e(V)
    matrix `Vw' = `Vw'[1..`k_count', 1..`k_count']
    matrix `chi2val' = `bw' * invsym(`Vw') * `bw''
    ereturn scalar chi2 = `chi2val'[1,1]
    ereturn scalar df_m = `k_count'
    if !missing(e(chi2)) {
        ereturn scalar p = chi2tail(`k_count', e(chi2))
    }
    else {
        ereturn scalar p = .
    }

    * --- Vuong test ---
    if "`vuong'" != "" & `ll_fenbreg' < . {
        ereturn scalar ll_fenbreg = `ll_fenbreg'
        local vuong_stat = sqrt(`N') * (`ll' - `ll_fenbreg') / `ll'
        * Proper Vuong stat computed in Mata and stored
    }

    * --- Display ---
    Display, level(`level') `irr' `vuong'

end

program define Display
    syntax [, Level(cilevel) IRR VUONG]

    local transform ""
    if "`irr'" != "" local transform "eform(IRR)"

    di _n as text e(title)
    di as text "Group variable: " as result e(ivar) ///
       as text _col(48) "Number of obs" _col(68) "= " as result %8.0fc e(N)
    di as text _col(48) "Number of groups" _col(68) "= " as result %8.0fc e(N_g)
    di _n as text "Log likelihood = " as result %12.4f e(ll) ///
       as text _col(48) "Wald chi2(" as result e(df_m) as text ")" ///
       _col(68) "= " as result %8.2f e(chi2)
    di as text "Dispersion:" _col(16) as result e(model) ///
       as text _col(48) "Prob > chi2" _col(68) "= " as result %8.4f e(p)
    di as text "ln(alpha) = " as result %9.6f e(lnalpha) ///
       as text "  alpha = " as result %9.6f e(alpha)
    if "`e(offset1)'" != "" {
        di as text "Offset: " as result "`e(offset1)'"
    }
    di as text "Inflate: " as result "`e(inflate)'"

    di ""
    di as text "{hline 13}{c TT}{hline 64}"
    di as text %12s "`e(depvar)'" " {c |}" _col(20) "Coefficient" ///
       _col(35) "Std. err." _col(49) "z" _col(55) "P>|z|" ///
       _col(65) "[`level'% conf. interval]"
    di as text "{hline 13}{c +}{hline 64}"

    _coef_table, level(`level') `transform'

    di as text "SEs clustered by " as result e(clustvar)

    * Vuong test
    if "`vuong'" != "" & e(ll_fenbreg) < . {
        di ""
        di as text "{hline 70}"
        di as text "Vuong test of ZINB vs NB (H0: models are equivalent)"
        di as text "  ZINB log-lik = " as result %12.4f e(ll) ///
           as text "    NB log-lik = " as result %12.4f e(ll_fenbreg)
        di as text "  Difference  = " as result %12.4f (e(ll) - e(ll_fenbreg))
        if e(ll) > e(ll_fenbreg) {
            di as text "  ZINB preferred (higher log-likelihood)"
        }
        else {
            di as text "  Standard NB preferred"
        }
        di as text "{hline 70}"
    }
end


// ============================================================
// Mata functions — embedded, compiled on first run
// ============================================================

mata:
mata set matastrict on

// --------------------------------------------------
// Main estimation driver
// --------------------------------------------------
void _zifenbreg_estimate(
    string scalar depname,
    string scalar indepnames,
    string scalar inflatenames,
    string scalar offname,
    string scalar ivarname,
    string scalar grpname,
    string scalar tousename,
    real scalar nb1,
    string scalar clustvarname,
    real scalar tol,
    real scalar maxiter,
    real scalar verbose,
    string scalar b0_count_name,
    string scalar b0_inflate_name,
    real scalar lnalpha_init,
    string scalar bname,
    string scalar biname,
    string scalar Vname,
    string scalar lnaname,
    string scalar llname,
    string scalar ll0name,
    string scalar convname,
    string scalar ll_fenbreg_name
)
{
    real colvector y, off, id, grp, clust
    real matrix X, Z
    real scalar N, k, k_z, N_g, i, i1, i2, ybar, denom

    st_view(y,   ., depname,      tousename)
    st_view(X,   ., indepnames,   tousename)
    st_view(off, ., offname,      tousename)
    st_view(id,  ., ivarname,     tousename)
    st_view(grp, ., grpname,      tousename)

    // Z matrix: inflate covariates + constant
    real matrix Zraw
    st_view(Zraw, ., inflatenames, tousename)
    Z = Zraw, J(rows(y), 1, 1)

    if (clustvarname == ivarname) {
        clust = grp
    } else {
        st_view(clust, ., clustvarname, tousename)
    }

    N = rows(y)
    k = cols(X)
    k_z = cols(Z)
    N_g = max(grp)

    real matrix info
    info = panelsetup(grp, 1)

    // --- Initialization ---
    real colvector beta, gamma, ai, yi, xbi
    real scalar lnalpha

    beta = st_matrix(b0_count_name)'
    if (rows(beta) > k) beta = beta[1::k]
    if (rows(beta) < k) beta = J(k, 1, 0)

    gamma = st_matrix(b0_inflate_name)'
    if (rows(gamma) > k_z) gamma = gamma[1::k_z]
    if (rows(gamma) < k_z) gamma = J(k_z, 1, 0)

    lnalpha = lnalpha_init

    // Initial unit effects
    ai = J(N_g, 1, 0)
    for (i = 1; i <= N_g; i++) {
        i1 = info[i, 1]
        i2 = info[i, 2]
        yi  = y[i1::i2]
        xbi = X[i1::i2, .] * beta + off[i1::i2]
        xbi = rowmin((xbi, J(rows(xbi), 1, 500)))
        ybar = mean(yi)
        denom = mean(exp(xbi))
        if (ybar > 0 & denom > 0 & denom < .) {
            ai[i] = ln(ybar / denom)
        }
    }

    // Initial inner solve
    _zifenbreg_solve_ci(y, X, Z, off, info, N_g, beta, gamma, lnalpha, ai, nb1, tol, 200)

    // --- Outer optimization ---
    real scalar ll_val, converged_flag
    real colvector theta

    theta = beta \ gamma \ lnalpha

    transmorphic M
    M = optimize_init()
    if (nb1) {
        optimize_init_evaluator(M, &_zifenbreg_conc_ll_d0())
        optimize_init_evaluatortype(M, "d0")
    }
    else {
        optimize_init_evaluator(M, &_zifenbreg_conc_ll())
        optimize_init_evaluatortype(M, "d1")
    }
    optimize_init_params(M, theta')
    optimize_init_argument(M, 1, y)
    optimize_init_argument(M, 2, X)
    optimize_init_argument(M, 3, Z)
    optimize_init_argument(M, 4, off)
    optimize_init_argument(M, 5, info)
    optimize_init_argument(M, 6, N_g)
    optimize_init_argument(M, 7, nb1)
    optimize_init_argument(M, 8, tol)
    optimize_init_argument(M, 9, &ai)
    optimize_init_technique(M, nb1 ? "bfgs" : "nr bfgs")
    optimize_init_singularHmethod(M, "hybrid")
    optimize_init_conv_ptol(M, tol)
    optimize_init_conv_vtol(M, tol)
    optimize_init_conv_nrtol(M, tol)
    optimize_init_tracelevel(M, verbose ? "value" : "none")
    optimize_init_conv_maxiter(M, maxiter)

    (void) _optimize(M)
    converged_flag = (optimize_result_errorcode(M) == 0)
    ll_val = optimize_result_value(M)

    theta = optimize_result_params(M)'
    beta = theta[1::k]
    gamma = theta[k+1::k+k_z]
    lnalpha = theta[k+k_z+1]

    // Final inner solve
    _zifenbreg_solve_ci(y, X, Z, off, info, N_g, beta, gamma, lnalpha, ai, nb1, tol, 200)

    // --- Null model log-likelihood (beta=0, gamma=0, alpha only) ---
    real scalar ll_0
    real colvector beta0, gamma0, ai0
    real scalar lna0
    beta0 = J(k, 1, 0)
    gamma0 = J(k_z, 1, 0)
    lna0 = lnalpha
    ai0 = ai
    _zifenbreg_solve_ci(y, X, Z, off, info, N_g, beta0, gamma0, lna0, ai0, nb1, tol, 200)
    ll_0 = _zifenbreg_ll_eval(y, X, Z, off, info, N_g, beta0, gamma0, lna0, ai0, nb1)

    // --- Sandwich variance ---
    real matrix Vmat, Hinv
    Hinv = optimize_result_V(M)
    Vmat = _zifenbreg_sandwich(y, X, Z, off, info, N_g, beta, gamma, lnalpha,
                                ai, nb1, clust, Hinv)

    // --- NB log-likelihood (no ZI) for Vuong test ---
    real scalar ll_fenbreg
    ll_fenbreg = _zifenbreg_nb_ll_eval(y, X, off, info, N_g, beta, lnalpha, ai, nb1)

    // --- Return to Stata ---
    st_matrix(bname, beta')
    st_matrix(biname, gamma')
    st_matrix(Vname, Vmat)
    st_numscalar(lnaname, lnalpha)
    st_numscalar(llname, ll_val)
    st_numscalar(ll0name, ll_0)
    st_numscalar(convname, converged_flag)
    st_numscalar(ll_fenbreg_name, ll_fenbreg)
}

// --------------------------------------------------
// Compute NB(0|mu,alpha) probability (log scale)
// Returns log P(y=0) under NB for a vector of mu values
// --------------------------------------------------
real colvector _zifenbreg_lnb0(
    real colvector mu,
    real scalar lnalpha,
    real scalar nb1
)
{
    real scalar alpha_val, kappa
    real colvector r, lp

    alpha_val = exp(lnalpha)

    if (!nb1) {
        // NB2: P(0) = (kappa/(kappa+mu))^kappa, kappa=1/alpha
        kappa = 1 / alpha_val
        lp = kappa :* (ln(kappa) :- ln(kappa :+ mu))
    }
    else {
        // NB1: P(0) = p^r, r=mu/alpha, p=1/(1+alpha)
        r = mu :/ alpha_val
        lp = r :* ln(1 / (1 + alpha_val))
    }
    return(lp)
}

// --------------------------------------------------
// Compute posterior weights for ZI model
// w_it = P(from count | y_it, params)
// tau_it = P(structural zero | y_it=0, params)
// --------------------------------------------------
void _zifenbreg_weights(
    real colvector y,
    real colvector mu,
    real colvector pi_vec,
    real scalar lnalpha,
    real scalar nb1,
    real colvector w,
    real colvector tau
)
{
    real scalar n
    real colvector lnb0, nb0, denom_vec

    n = rows(y)
    w = J(n, 1, 1)
    tau = J(n, 1, 0)

    lnb0 = _zifenbreg_lnb0(mu, lnalpha, nb1)
    nb0 = exp(lnb0)

    real scalar j
    for (j = 1; j <= n; j++) {
        if (y[j] == 0) {
            denom_vec = pi_vec[j] + (1 - pi_vec[j]) * nb0[j]
            if (denom_vec > 1e-300) {
                tau[j] = pi_vec[j] / denom_vec
                w[j] = (1 - pi_vec[j]) * nb0[j] / denom_vec
            }
            else {
                tau[j] = 0.5
                w[j] = 0.5
            }
        }
    }
}

// --------------------------------------------------
// Inner NR loop: solve for unit effects given (beta, gamma, lnalpha)
// Modified for ZI weights on y=0 observations
// --------------------------------------------------
void _zifenbreg_solve_ci(
    real colvector y,
    real matrix X,
    real matrix Z,
    real colvector off,
    real matrix info,
    real scalar N_g,
    real colvector beta,
    real colvector gamma,
    real scalar lnalpha,
    real colvector ai,
    real scalar nb1,
    real scalar tol,
    real scalar maxiter
)
{
    real scalar i, i1, i2, iter, alpha_val, kappa, score_i, hess_i, step, pp
    real colvector yi, mu_i, xbi, r_i, pi_i, w_i, tau_i, zgi

    alpha_val = exp(lnalpha)

    for (i = 1; i <= N_g; i++) {
        i1 = info[i, 1]
        i2 = info[i, 2]
        yi  = y[i1::i2]
        xbi = X[i1::i2, .] * beta + off[i1::i2]
        zgi = Z[i1::i2, .] * gamma
        pi_i = invlogit(zgi)

        for (iter = 1; iter <= maxiter; iter++) {
            mu_i = exp(rowmin((ai[i] :+ xbi, J(rows(xbi), 1, 500))))

            // Compute ZI weights
            _zifenbreg_weights(yi, mu_i, pi_i, lnalpha, nb1, w_i, tau_i)

            if (!nb1) {
                kappa = 1 / alpha_val
                score_i = sum(w_i :* kappa :* (yi :- mu_i) :/ (kappa :+ mu_i))
                hess_i  = -sum(w_i :* kappa :* (yi :+ kappa) :* mu_i :/
                           (kappa :+ mu_i):^2)
            }
            else {
                r_i = mu_i :/ alpha_val
                pp = 1 / (1 + alpha_val)
                score_i = sum(w_i :* r_i :* (digamma(yi :+ r_i) :- digamma(r_i) :+ ln(pp)))
                // Use only the guaranteed-negative trigamma term as Hessian
                // (the full Hessian = score + this term, but score can make it positive)
                hess_i = sum(w_i :* r_i:^2 :*
                    (trigamma(yi :+ r_i) :- trigamma(r_i)))
            }

            if (hess_i >= 0) hess_i = -abs(score_i) - 1
            if (score_i >= . | hess_i >= .) break
            step = -score_i / hess_i
            if (step > 10)  step = 10
            if (step < -10) step = -10
            ai[i] = ai[i] + step

            if (abs(step) < tol) break
        }
    }
}

// --------------------------------------------------
// Evaluate full ZI NB log-likelihood
// --------------------------------------------------
real scalar _zifenbreg_ll_eval(
    real colvector y,
    real matrix X,
    real matrix Z,
    real colvector off,
    real matrix info,
    real scalar N_g,
    real colvector beta,
    real colvector gamma,
    real scalar lnalpha,
    real colvector ai,
    real scalar nb1
)
{
    real scalar ll, alpha_val, kappa, i, i1, i2, j
    real colvector yi, mu_i, xbi, r_i, pi_i, zgi, lnb0_i, ll_count_i

    alpha_val = exp(lnalpha)
    ll = 0

    for (i = 1; i <= N_g; i++) {
        i1 = info[i, 1]
        i2 = info[i, 2]
        yi  = y[i1::i2]
        xbi = X[i1::i2, .] * beta + off[i1::i2]
        mu_i = exp(rowmin((ai[i] :+ xbi, J(rows(xbi), 1, 500))))
        zgi = Z[i1::i2, .] * gamma
        pi_i = invlogit(zgi)

        if (!nb1) {
            kappa = max((1 / alpha_val, 1e-10))
            ll_count_i = lngamma(yi :+ kappa) :- lngamma(kappa) :- lngamma(yi :+ 1) :+
                kappa :* ln(kappa) :- kappa :* ln(kappa :+ mu_i) :+
                yi :* ln(mu_i :+ 1e-300) :- yi :* ln(kappa :+ mu_i)
        }
        else {
            r_i = rowmax((mu_i :/ alpha_val, J(rows(mu_i), 1, 1e-10)))
            ll_count_i = lngamma(yi :+ r_i) :- lngamma(r_i) :- lngamma(yi :+ 1) :+
                r_i :* ln(1 / (1 + alpha_val)) :+
                yi :* ln(alpha_val / (1 + alpha_val))
        }

        // ZI log-likelihood
        for (j = 1; j <= rows(yi); j++) {
            if (yi[j] == 0) {
                // log(pi + (1-pi)*NB(0))
                ll = ll + ln(pi_i[j] + (1 - pi_i[j]) * exp(ll_count_i[j]) + 1e-300)
            }
            else {
                // log((1-pi)*NB(y))
                ll = ll + ln(1 - pi_i[j] + 1e-300) + ll_count_i[j]
            }
        }
    }

    if (ll >= .) ll = -1e100
    return(ll)
}

// --------------------------------------------------
// d0 evaluator: concentrated ZI NB log-likelihood only (for NB1)
// --------------------------------------------------
void _zifenbreg_conc_ll_d0(
    real scalar todo,
    real rowvector params,
    real colvector y,
    real matrix X,
    real matrix Z,
    real colvector off,
    real matrix info,
    real scalar N_g,
    real scalar nb1,
    real scalar tol,
    pointer(real colvector) scalar pai,
    real scalar ll,
    real rowvector g,
    real matrix H
)
{
    real scalar k, k_z, lnalpha
    real colvector beta, gamma, ai

    k = cols(X)
    k_z = cols(Z)
    beta = params[1::k]'
    gamma = params[k+1::k+k_z]'
    lnalpha = params[k+k_z+1]
    if (lnalpha < -10) lnalpha = -10
    if (lnalpha > 10)  lnalpha = 10

    ai = *pai
    _zifenbreg_solve_ci(y, X, Z, off, info, N_g, beta, gamma, lnalpha, ai, nb1, tol, 200)
    *pai = ai

    ll = _zifenbreg_ll_eval(y, X, Z, off, info, N_g, beta, gamma, lnalpha, ai, nb1)

    if (ll >= . | ll > 0 | ll == 0) {
        ll = -1e100
    }
}

// --------------------------------------------------
// d1 evaluator: concentrated ZI NB log-likelihood + gradient
// --------------------------------------------------
void _zifenbreg_conc_ll(
    real scalar todo,
    real rowvector params,
    real colvector y,
    real matrix X,
    real matrix Z,
    real colvector off,
    real matrix info,
    real scalar N_g,
    real scalar nb1,
    real scalar tol,
    pointer(real colvector) scalar pai,
    real scalar ll,
    real rowvector g,
    real matrix H
)
{
    real scalar k, k_z, alpha_val, kappa, lnalpha, pp, i, i1, i2
    real colvector beta, gamma, ai, gbeta, ggamma, yi, mu_i, xbi, zgi
    real colvector w_i, tau_i, pi_i, r_i, psi_diff, dlda
    real scalar glnalpha

    k = cols(X)
    k_z = cols(Z)
    beta = params[1::k]'
    gamma = params[k+1::k+k_z]'
    lnalpha = params[k+k_z+1]
    if (lnalpha < -10) lnalpha = -10
    if (lnalpha > 10)  lnalpha = 10
    alpha_val = exp(lnalpha)

    ai = *pai

    // Solve inner problem
    _zifenbreg_solve_ci(y, X, Z, off, info, N_g, beta, gamma, lnalpha, ai, nb1, tol, 200)
    *pai = ai

    // Compute log-likelihood
    ll = _zifenbreg_ll_eval(y, X, Z, off, info, N_g, beta, gamma, lnalpha, ai, nb1)

    if (ll >= . | ll > 0 | ll == 0) {
        ll = -1e100
        if (todo >= 1) g = J(1, k + k_z + 1, 0)
        return
    }

    if (todo == 0) return

    // Gradient via envelope theorem
    gbeta = J(k, 1, 0)
    ggamma = J(k_z, 1, 0)
    glnalpha = 0

    for (i = 1; i <= N_g; i++) {
        i1 = info[i, 1]
        i2 = info[i, 2]
        yi  = y[i1::i2]
        xbi = X[i1::i2, .] * beta + off[i1::i2]
        mu_i = exp(rowmin((ai[i] :+ xbi, J(rows(xbi), 1, 500))))
        zgi = Z[i1::i2, .] * gamma
        pi_i = invlogit(zgi)

        // Compute posterior weights
        _zifenbreg_weights(yi, mu_i, pi_i, lnalpha, nb1, w_i, tau_i)

        // --- d(LL)/d(beta) --- weighted by w_it for y=0 terms
        if (!nb1) {
            kappa = 1 / alpha_val
            gbeta = gbeta + X[i1::i2, .]' * (w_i :* kappa :* (yi :- mu_i) :/ (kappa :+ mu_i))

            // --- d(LL)/d(lnalpha) --- weighted version
            glnalpha = glnalpha - kappa * sum(w_i :* (
                digamma(yi :+ kappa) :- digamma(kappa) :+
                ln(kappa) :+ 1 :- ln(kappa :+ mu_i) :-
                (yi :+ kappa) :/ (kappa :+ mu_i)
            ))
        }
        else {
            r_i = mu_i :/ alpha_val
            pp = 1 / (1 + alpha_val)
            psi_diff = digamma(yi :+ r_i) :- digamma(r_i)

            gbeta = gbeta + X[i1::i2, .]' * (w_i :* (1 / alpha_val) :* (psi_diff :+ ln(pp)) :* mu_i)

            dlda = (-mu_i :/ alpha_val^2) :* (psi_diff :+ ln(pp)) :+
                   r_i :* (-1 / (1 + alpha_val)) :+
                   yi :* (1 / (1 + alpha_val))
            glnalpha = glnalpha + alpha_val * sum(w_i :* dlda)
        }

        // --- d(LL)/d(gamma) ---
        // y=0: d/dgamma log(pi+(1-pi)*NB0) = (1-NB0)*pi*(1-pi) / f * z
        // y>0: d/dgamma log(1-pi) = -pi * z
        real colvector g_inflate_i
        real scalar j
        g_inflate_i = J(rows(yi), 1, 0)

        real colvector lnb0_i, nb0_i, f_i
        lnb0_i = _zifenbreg_lnb0(mu_i, lnalpha, nb1)
        nb0_i = exp(lnb0_i)

        for (j = 1; j <= rows(yi); j++) {
            if (yi[j] == 0) {
                f_i = pi_i[j] + (1 - pi_i[j]) * nb0_i[j]
                if (f_i > 1e-300) {
                    g_inflate_i[j] = (1 - nb0_i[j]) * pi_i[j] * (1 - pi_i[j]) / f_i
                }
            }
            else {
                g_inflate_i[j] = -pi_i[j]
            }
        }

        ggamma = ggamma + Z[i1::i2, .]' * g_inflate_i
    }

    g = (gbeta \ ggamma \ glnalpha)'
}

// --------------------------------------------------
// Sandwich variance estimator (panel-clustered)
// --------------------------------------------------
real matrix _zifenbreg_sandwich(
    real colvector y,
    real matrix X,
    real matrix Z,
    real colvector off,
    real matrix info,
    real scalar N_g,
    real colvector beta,
    real colvector gamma,
    real scalar lnalpha,
    real colvector ai,
    real scalar nb1,
    real colvector clust,
    real matrix Hinv
)
{
    real scalar k, k_z, alpha_val, kappa, pp, N, i, i1, i2, nc, j, jj
    real colvector yi, mu_i, xbi, zgi, pi_i, w_i, tau_i, ga_i, r_i, psi_diff, dlda
    real colvector uclust, sel, lnb0_i, nb0_i
    real matrix scores, S, B, V
    real scalar f_ij

    k = cols(X)
    k_z = cols(Z)
    alpha_val = exp(lnalpha)
    N = rows(y)

    // Per-observation score contributions: N x (k + k_z + 1)
    scores = J(N, k + k_z + 1, 0)

    for (i = 1; i <= N_g; i++) {
        i1 = info[i, 1]
        i2 = info[i, 2]
        yi  = y[i1::i2]
        xbi = X[i1::i2, .] * beta + off[i1::i2]
        mu_i = exp(rowmin((ai[i] :+ xbi, J(rows(xbi), 1, 500))))
        zgi = Z[i1::i2, .] * gamma
        pi_i = invlogit(zgi)

        _zifenbreg_weights(yi, mu_i, pi_i, lnalpha, nb1, w_i, tau_i)

        lnb0_i = _zifenbreg_lnb0(mu_i, lnalpha, nb1)
        nb0_i = exp(lnb0_i)

        if (!nb1) {
            kappa = 1 / alpha_val
            // Beta scores
            scores[i1::i2, 1::k] = X[i1::i2, .] :* (w_i :* kappa :* (yi :- mu_i) :/ (kappa :+ mu_i))

            // lnalpha scores
            ga_i = -kappa :* w_i :* (
                digamma(yi :+ kappa) :- digamma(kappa) :+
                ln(kappa) :+ 1 :- ln(kappa :+ mu_i) :-
                (yi :+ kappa) :/ (kappa :+ mu_i)
            )
            scores[i1::i2, k+k_z+1] = ga_i
        }
        else {
            r_i = mu_i :/ alpha_val
            pp = 1 / (1 + alpha_val)
            psi_diff = digamma(yi :+ r_i) :- digamma(r_i)

            scores[i1::i2, 1::k] = X[i1::i2, .] :*
                (w_i :* (1/alpha_val) :* (psi_diff :+ ln(pp)) :* mu_i)

            dlda = (-mu_i :/ alpha_val^2) :* (psi_diff :+ ln(pp)) :+
                   r_i :* (-1 / (1 + alpha_val)) :+
                   yi :* (1 / (1 + alpha_val))
            scores[i1::i2, k+k_z+1] = alpha_val :* w_i :* dlda
        }

        // Inflate equation scores
        real colvector g_inf_i
        g_inf_i = J(rows(yi), 1, 0)
        for (jj = 1; jj <= rows(yi); jj++) {
            if (yi[jj] == 0) {
                f_ij = pi_i[jj] + (1 - pi_i[jj]) * nb0_i[jj]
                if (f_ij > 1e-300) {
                    g_inf_i[jj] = (1 - nb0_i[jj]) * pi_i[jj] * (1 - pi_i[jj]) / f_ij
                }
            }
            else {
                g_inf_i[jj] = -pi_i[jj]
            }
        }
        scores[i1::i2, k+1::k+k_z] = Z[i1::i2, .] :* g_inf_i
    }

    // Cluster-level score sums
    uclust = uniqrows(clust)
    nc = rows(uclust)
    S = J(nc, k + k_z + 1, 0)

    for (j = 1; j <= nc; j++) {
        sel = selectindex(clust :== uclust[j])
        S[j, .] = colsum(scores[sel, .])
    }

    // Meat and sandwich
    B = S' * S
    V = Hinv * B * Hinv

    // Small-sample correction
    real scalar ssc
    ssc = nc / (nc - 1)
    V = ssc :* V

    return(V)
}

// --------------------------------------------------
// Standard NB log-likelihood (no ZI) for Vuong test
// --------------------------------------------------
real scalar _zifenbreg_nb_ll_eval(
    real colvector y,
    real matrix X,
    real colvector off,
    real matrix info,
    real scalar N_g,
    real colvector beta,
    real scalar lnalpha,
    real colvector ai,
    real scalar nb1
)
{
    real scalar ll, alpha_val, kappa, i, i1, i2
    real colvector yi, mu_i, xbi, r_i

    alpha_val = exp(lnalpha)
    ll = 0

    for (i = 1; i <= N_g; i++) {
        i1 = info[i, 1]
        i2 = info[i, 2]
        yi  = y[i1::i2]
        xbi = X[i1::i2, .] * beta + off[i1::i2]
        mu_i = exp(rowmin((ai[i] :+ xbi, J(rows(xbi), 1, 500))))

        if (!nb1) {
            kappa = max((1 / alpha_val, 1e-10))
            ll = ll + sum(
                lngamma(yi :+ kappa) :- lngamma(kappa) :- lngamma(yi :+ 1) :+
                kappa :* ln(kappa) :- kappa :* ln(kappa :+ mu_i) :+
                yi :* ln(mu_i :+ 1e-300) :- yi :* ln(kappa :+ mu_i)
            )
        }
        else {
            r_i = rowmax((mu_i :/ alpha_val, J(rows(mu_i), 1, 1e-10)))
            ll = ll + sum(
                lngamma(yi :+ r_i) :- lngamma(r_i) :- lngamma(yi :+ 1) :+
                r_i :* ln(1 / (1 + alpha_val)) :+
                yi :* ln(alpha_val / (1 + alpha_val))
            )
        }
    }
    if (ll >= .) ll = -1e100
    return(ll)
}

end
