*! ebdid v1.0.0
*! Empirical Bayes Shrinkage Difference-in-Differences

program define ebdid, eclass sortpreserve
    version 15.0

    syntax varname(numeric) [if] [in], ///
        Group(varname numeric) Time(varname numeric) ///
        TReated(varname numeric) ///
        [Breps(integer 200) VCalib(string) Level(cilevel)]

    * --- defaults ---
    if "`vcalib'" == "" local vcalib "untreated"
    if !inlist("`vcalib'", "untreated", "treated", "both") {
        di as error "vcalib() must be: untreated, treated, or both"
        exit 198
    }
    if "`level'" == "" local level = c(level)
    local z = invnormal(1 - (1 - `level'/100)/2)

    marksample touse
    local depvar `varlist'
    local grp_label : value label `group'

    preserve
    quietly {
        keep if `touse'

        * ── 1. IDENTIFY GROUPS ──
        tempvar ever_tr tstart gtag
        bysort `group': egen byte `ever_tr' = max(`treated')
        bysort `group': egen double `tstart' = ///
            min(cond(`treated' == 1, `time', .))
        egen `gtag' = tag(`group')

        count if `gtag' & `ever_tr'
        local G_treat = r(N)
        count if `gtag' & !`ever_tr'
        local G_ctrl = r(N)

        if `G_treat' < 4 {
            restore
            di as error "ebdid requires at least 4 treated groups"
            exit 198
        }
        if `G_ctrl' < 1 {
            restore
            di as error ///
                "ebdid requires at least 1 never-treated group"
            exit 198
        }

        su `time', meanonly
        local T_min = r(min)
        local T_max = r(max)

        * ── 2. COLLAPSE TO GROUP × TIME ──
        collapse (mean) `depvar' (first) `ever_tr' `tstart', ///
            by(`group' `time')

        * store treated group IDs and treatment starts
        levelsof `group' if `ever_tr', local(tr_groups)
        local gi = 0
        foreach g of local tr_groups {
            local ++gi
            su `tstart' if `group' == `g', meanonly
            local tstart_`gi' = r(mean)
            local gid_`gi' = `g'
        }

        * ── 3. PLACEBOS FOR V ──
        sort `group' `time'
        by `group': gen double _dy = `depvar' - `depvar'[_n-1]

        if "`vcalib'" != "treated" {
            * untreated LOO placebos (all periods)
            bysort `time': egen double _sdy_c = ///
                total(cond(`ever_tr' == 0, _dy, 0))
            bysort `time': egen double _nc = ///
                total(`ever_tr' == 0 & _dy < .)
            gen double _dy_loo = (_sdy_c - _dy) / (_nc - 1) ///
                if `ever_tr' == 0 & _nc > 1
            gen double _pc = _dy - _dy_loo ///
                if `ever_tr' == 0 & _dy < . & _dy_loo < .
            gen double _pc2 = _pc * _pc
            su _pc2, meanonly
            local msp_c = r(mean)
        }

        if "`vcalib'" != "untreated" {
            * treated pre-period placebos
            bysort `time': egen double _dy_cm = ///
                mean(cond(`ever_tr' == 0, _dy, .))
            gen double _pt = _dy - _dy_cm ///
                if `ever_tr' == 1 & _dy < . & `time' < `tstart'
            gen double _pt2 = _pt * _pt
            su _pt2, meanonly
            local msp_t = r(mean)
        }

        * V_scale per treated group
        local sum_vs = 0
        forvalues i = 1/`G_treat' {
            local Tpre = `tstart_`i'' - `T_min'
            local Tpost = `T_max' - `tstart_`i'' + 1
            if `Tpre' < 1 | `Tpost' < 1 {
                restore
                di as error ///
                    "group `gid_`i'': no pre or post periods"
                exit 198
            }
            local vs_`i' = (1/`Tpost' + 1/`Tpre') / 2
            local sum_vs = `sum_vs' + `vs_`i''
        }
        local mean_vs = `sum_vs' / `G_treat'

        if "`vcalib'" == "untreated" {
            local V_hat = `msp_c' * `mean_vs'
        }
        else if "`vcalib'" == "treated" {
            local V_hat = `msp_t' * `mean_vs'
        }
        else {
            su _pc2, meanonly
            local wc = r(N)
            su _pt2, meanonly
            local wt = r(N)
            local V_hat = ///
                (`msp_c'*`wc' + `msp_t'*`wt') ///
                / (`wc' + `wt') * `mean_vs'
        }

        * ── 4. GROUP-LEVEL DiDs ──
        matrix _ebdid_tols = J(`G_treat', 1, .)
        matrix _ebdid_gid  = J(`G_treat', 1, .)

        forvalues i = 1/`G_treat' {
            local g = `gid_`i''
            local ts = `tstart_`i''

            su `depvar' if `group'==`g' & `time'>=`ts', meanonly
            local yg1 = r(mean)
            su `depvar' if `group'==`g' & `time'<`ts', meanonly
            local yg0 = r(mean)
            su `depvar' if `ever_tr'==0 & `time'>=`ts', meanonly
            local yc1 = r(mean)
            su `depvar' if `ever_tr'==0 & `time'<`ts', meanonly
            local yc0 = r(mean)

            matrix _ebdid_tols[`i',1] = (`yg1'-`yg0')-(`yc1'-`yc0')
            matrix _ebdid_gid[`i',1]  = `g'
        }

        * ── 5. JS SHRINKAGE ──
        local sum_tau = 0
        forvalues i = 1/`G_treat' {
            local sum_tau = `sum_tau' + _ebdid_tols[`i',1]
        }
        local tau_bar = `sum_tau' / `G_treat'

        local S2 = 0
        forvalues i = 1/`G_treat' {
            local S2 = `S2' + (_ebdid_tols[`i',1] - `tau_bar')^2
        }

        local B_hat = 0
        if `S2' > 0 & `G_treat' >= 4 {
            local B_hat = max(0, ///
                1 - (`G_treat' - 3) * `V_hat' / `S2')
        }

        matrix _ebdid_teb = J(`G_treat', 1, .)
        forvalues i = 1/`G_treat' {
            matrix _ebdid_teb[`i',1] = `tau_bar' + ///
                `B_hat' * (_ebdid_tols[`i',1] - `tau_bar')
        }

        * ── 6. BOOTSTRAP ──
        tempname sc_V sc_B
        scalar `sc_V' = `V_hat'
        scalar `sc_B' = `breps'
        mata: st_matrix("_ebdid_se", ///
            _ebdid_boot_se( ///
                st_matrix("_ebdid_tols"), ///
                st_numscalar("`sc_V'"), ///
                st_numscalar("`sc_B'")))

        * ── BUILD e(b) AND e(V) ──
        matrix _ebdid_b = _ebdid_teb'
        matrix _ebdid_V = J(`G_treat', `G_treat', 0)
        forvalues i = 1/`G_treat' {
            matrix _ebdid_V[`i',`i'] = _ebdid_se[`i',1]^2
        }

        local cnames ""
        forvalues i = 1/`G_treat' {
            local g = int(_ebdid_gid[`i',1])
            local cnames "`cnames' g`g':tau"
        }
        matrix colnames _ebdid_b = `cnames'
        matrix colnames _ebdid_V = `cnames'
        matrix rownames _ebdid_V = `cnames'
    }
    restore

    * ── POST RESULTS ──
    ereturn post _ebdid_b _ebdid_V
    ereturn local cmd "ebdid"
    ereturn local cmdline `"ebdid `0'"'
    ereturn local depvar "`depvar'"
    ereturn local vcalib "`vcalib'"
    ereturn local properties "b V"

    ereturn scalar B_hat   = `B_hat'
    ereturn scalar V_hat   = `V_hat'
    ereturn scalar tau_bar = `tau_bar'
    ereturn scalar S2      = `S2'
    ereturn scalar G_treat = `G_treat'
    ereturn scalar G_ctrl  = `G_ctrl'
    ereturn scalar breps   = `breps'
    ereturn scalar level   = `level'

    ereturn matrix tau_ols  = _ebdid_tols
    ereturn matrix group_id = _ebdid_gid
    ereturn matrix se_boot  = _ebdid_se

    * ── DISPLAY ──
    _ebdid_display, grplabel(`grp_label')
end


* ── DISPLAY SUBROUTINE ──
program _ebdid_display
    syntax [, grplabel(string)]

    local level = e(level)
    local z = invnormal(1 - (1 - `level'/100)/2)
    local G = e(G_treat)

    di ""
    di as txt "Empirical Bayes Shrinkage DiD" ///
        _col(49) "Number of groups"
    di as txt "{hline 48}" ///
        _col(49) "  Treated  = " as res %5.0f e(G_treat)
    di as txt "Shrinkage factor (B)  = " ///
        as res %9.4f e(B_hat) ///
        _col(49) as txt "  Control  = " as res %5.0f e(G_ctrl)
    di as txt "Noise variance (V)    = " ///
        as res %9.4f e(V_hat)
    di as txt "Grand mean (tau_bar)  = " ///
        as res %9.4f e(tau_bar) ///
        _col(49) as txt "  Bootstrap = " ///
        as res %5.0f e(breps) as txt " reps"
    di as txt "V calibration         = " ///
        as res "`e(vcalib)'"
    di ""

    di as txt "{hline 13}{c TT}{hline 64}"
    di as txt %12s "Group" " {c |}" ///
        %10s "tau_ols" %10s "tau_eb" ///
        %10s "Std.Err." %8s "z" %8s "P>|z|" ///
        "  [`level'% Conf. Interval]"
    di as txt "{hline 13}{c +}{hline 64}"

    tempname gid tols se
    matrix `gid'  = e(group_id)
    matrix `tols' = e(tau_ols)
    matrix `se'   = e(se_boot)

    forvalues i = 1/`G' {
        local g = int(`gid'[`i',1])

        * use value label if available
        if "`grplabel'" != "" {
            local gdisp : label `grplabel' `g'
        }
        else {
            local gdisp "`g'"
        }

        local t_ols = `tols'[`i',1]
        local t_eb  = el(e(b), 1, `i')
        local s     = `se'[`i',1]

        if `s' > 1e-10 {
            local zv = `t_eb' / `s'
            local pv = 2 * normal(-abs(`zv'))
        }
        else {
            local zv = 0
            local pv = 1
        }

        di as txt %12s "`gdisp'" " {c |}" ///
            as res %10.4f `t_ols' ///
            %10.4f `t_eb' ///
            %10.4f `s' ///
            %8.2f `zv' ///
            %8.3f `pv' ///
            %10.4f (`t_eb' - `z' * `s') ///
            %10.4f (`t_eb' + `z' * `s')
    }

    di as txt "{hline 13}{c BT}{hline 64}"
end


* ── MATA: PARAMETRIC BOOTSTRAP ──
mata:
real colvector _ebdid_boot_se(real colvector tau_hat, ///
    real scalar V, real scalar B_boot)
{
    G = rows(tau_hat)
    sdV = sqrt(V)
    tau_tilde = J(G, B_boot, .)

    for (b = 1; b <= B_boot; b++) {
        tau_star = tau_hat + rnormal(G, 1, 0, sdV)
        tbar_s = mean(tau_star)
        S2_s = sum((tau_star :- tbar_s):^2)

        B_s = 0
        if (S2_s > 0 & G >= 4) {
            B_s = max((0, 1 - (G - 3) * V / S2_s))
        }

        tau_tilde[., b] = tbar_s :+ B_s * (tau_star :- tbar_s)
    }

    se = J(G, 1, .)
    for (g = 1; g <= G; g++) {
        se[g] = sqrt(variance(tau_tilde[g, .]'))
    }
    return(se)
}
end
