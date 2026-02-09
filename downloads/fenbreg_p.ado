*! fenbreg_p v1.0.0 — Predict after fenbreg

program define fenbreg_p
    version 16.0

    syntax newvarname [if] [in] , [ MU XB FE N ]

    local type "`mu'`xb'`fe'`n'"
    if "`type'" == "" local type "mu"

    local nopt : word count `mu' `xb' `fe' `n'
    if `nopt' > 1 {
        di as error "only one of mu, xb, fe, or n allowed"
        exit 198
    }

    if "`e(cmd)'" != "fenbreg" {
        error 301
    }

    marksample touse, novarlist
    local depvar "`e(depvar)'"
    local ivar "`e(ivar)'"

    tempvar xbvar

    * Compute xb = X*beta (no constant)
    qui _predict double `xbvar' if `touse', xb

    if "`type'" == "xb" {
        gen `typlist' `varlist' = `xbvar' if `touse'
        label variable `varlist' "Linear prediction (xb)"
        exit
    }

    * Compute unit effects: a_i = log(c_i)
    * c_i = ybar_i / mean_t(exp(xb_it))
    tempvar fe_var mu_var grp expxb ybar denom
    sort `ivar'

    if "`e(offset1)'" != "" {
        * Account for offset
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

    if "`type'" == "n" | "`type'" == "mu" {
        if "`e(offset1)'" != "" {
            tempvar offv2
            if "`e(exposure)'" != "" {
                qui gen double `offv2' = ln(`e(exposure)') if `touse'
            }
            else if "`e(offset)'" != "" {
                qui gen double `offv2' = `e(offset)' if `touse'
            }
            else {
                qui gen double `offv2' = 0 if `touse'
            }
            qui gen double `mu_var' = exp(`fe_var' + `xbvar' + `offv2') if `touse'
        }
        else {
            qui gen double `mu_var' = exp(`fe_var' + `xbvar') if `touse'
        }

        gen `typlist' `varlist' = `mu_var' if `touse'
        if "`type'" == "mu" {
            label variable `varlist' "Predicted mean (mu_it)"
        }
        else {
            label variable `varlist' "Predicted count (n)"
        }
        exit
    }
end
