*! zifenbreg_p v1.0.0 — Predict after zifenbreg

program define zifenbreg_p
    version 16.0

    syntax newvarname [if] [in] , [ MU PR XB ZG FE YHAT ]

    local type "`mu'`pr'`xb'`zg'`fe'`yhat'"
    if "`type'" == "" local type "yhat"

    local nopt : word count `mu' `pr' `xb' `zg' `fe' `yhat'
    if `nopt' > 1 {
        di as error "only one of mu, pr, xb, zg, fe, or yhat allowed"
        exit 198
    }

    if "`e(cmd)'" != "zifenbreg" {
        error 301
    }

    marksample touse, novarlist
    local depvar "`e(depvar)'"
    local ivar "`e(ivar)'"
    local inflate "`e(inflate)'"
    local k_count = e(k_count)
    local k_inflate = e(k_inflate)

    * --- Count equation: xb = X*beta (no constant) ---
    if "`type'" == "xb" {
        tempvar xbvar
        * Get count betas from e(b)
        tempname b_count
        matrix `b_count' = e(b)
        matrix `b_count' = `b_count'[1, 1..`k_count']
        qui matrix score double `xbvar' = `b_count' if `touse'
        gen `typlist' `varlist' = `xbvar' if `touse'
        label variable `varlist' "Linear prediction, count equation (xb)"
        exit
    }

    * --- Inflate equation: zg = Z*gamma (with constant) ---
    if "`type'" == "zg" {
        tempvar zgvar
        tempname b_inflate
        matrix `b_inflate' = e(b)
        matrix `b_inflate' = `b_inflate'[1, `k_count'+1..`k_count'+`k_inflate']
        qui matrix score double `zgvar' = `b_inflate' if `touse'
        gen `typlist' `varlist' = `zgvar' if `touse'
        label variable `varlist' "Linear prediction, inflate equation (zg)"
        exit
    }

    * --- P(inflate) = invlogit(zg) ---
    if "`type'" == "pr" {
        tempvar zgvar
        tempname b_inflate
        matrix `b_inflate' = e(b)
        matrix `b_inflate' = `b_inflate'[1, `k_count'+1..`k_count'+`k_inflate']
        qui matrix score double `zgvar' = `b_inflate' if `touse'
        gen `typlist' `varlist' = invlogit(`zgvar') if `touse'
        label variable `varlist' "P(inflate), logit prediction"
        exit
    }

    * --- Compute xb for remaining options ---
    tempvar xbvar
    tempname b_count
    matrix `b_count' = e(b)
    matrix `b_count' = `b_count'[1, 1..`k_count']
    qui matrix score double `xbvar' = `b_count' if `touse'

    * --- Unit fixed effect ---
    tempvar fe_var expxb ybar denom
    sort `ivar'

    if "`e(offset1)'" != "" {
        tempvar offv
        if "`e(exposure)'" != "" {
            qui gen double `offv' = ln(`e(exposure)') if `touse'
        }
        else if "`e(offset)'" != "" {
            qui gen double `offv' = `e(offset)' if `touse'
        }
        else {
            qui gen double `offv' = 0 if `touse'
        }
        qui gen double `expxb' = exp(`xbvar' + `offv') if `touse'
    }
    else {
        qui gen double `expxb' = exp(`xbvar') if `touse'
    }

    qui by `ivar': egen double `ybar' = mean(`depvar') if `touse'
    qui by `ivar': egen double `denom' = mean(`expxb') if `touse'
    qui gen double `fe_var' = ln(`ybar' / `denom') if `touse' & `ybar' > 0
    qui replace `fe_var' = 0 if `touse' & `ybar' == 0

    if "`type'" == "fe" {
        gen `typlist' `varlist' = `fe_var' if `touse'
        label variable `varlist' "Unit fixed effect (log c_i)"
        exit
    }

    * --- mu = exp(fe + xb + offset) ---
    tempvar mu_var
    if "`e(offset1)'" != "" {
        qui gen double `mu_var' = exp(`fe_var' + `xbvar' + `offv') if `touse'
    }
    else {
        qui gen double `mu_var' = exp(`fe_var' + `xbvar') if `touse'
    }

    if "`type'" == "mu" {
        gen `typlist' `varlist' = `mu_var' if `touse'
        label variable `varlist' "Predicted count mean (mu_it)"
        exit
    }

    * --- yhat = (1-pi)*mu ---
    if "`type'" == "yhat" {
        tempvar zgvar pi_var
        tempname b_inflate
        matrix `b_inflate' = e(b)
        matrix `b_inflate' = `b_inflate'[1, `k_count'+1..`k_count'+`k_inflate']
        qui matrix score double `zgvar' = `b_inflate' if `touse'
        qui gen double `pi_var' = invlogit(`zgvar') if `touse'
        gen `typlist' `varlist' = (1 - `pi_var') * `mu_var' if `touse'
        label variable `varlist' "Predicted mean E[y] = (1-pi)*mu"
        exit
    }
end
