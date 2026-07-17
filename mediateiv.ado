*! mediateiv.ado v3.0  16jul2026
*! IV causal mediation analysis with D*M interaction.
*! Supports both endogenous and exogenous treatment.
*!
*! v3.0 changes:
*!   - New ivregress-style syntax:
*!       mediateiv y (d=z1) (m=z2) x1 x2 x3        (D endogenous, 3-stage)
*!       mediateiv y d (m=z3) x1 x2 x3              (D exogenous, 2-stage)
*!   - When D exogenous: 2-stage estimation (no D instrumentation)
*!   - Factor variable notation in controls (i.prov)
*!   - Bootstrap SE/CI/p-value, N at table bottom
*!
*! Original author: Xiliang Zhao (zhaoxiliang@gmail.com), 2024

program mediateiv, rclass
    version 19
    syntax anything(equalok everything) [if] [in], [REPS(integer 100) SEED(string)]

    *========================================================
    * Parse: y (d=ivd) (m=ivm) controls
    *     or y d (m=ivm) controls
    *========================================================
    gettoken yvar rest : anything
    confirm variable `yvar'

    local d_endog  = 0
    local dvar     ""
    local ivd      ""
    local mvar     ""
    local ivm      ""
    local controls ""

    while `"`rest'"' != "" {
        gettoken tok rest : rest, bind
        if regexm(`"`tok'"', "^\((.*)\)$") {
            local content = regexs(1)
            local eqpos = strpos(`"`content'"', "=")
            if `eqpos' > 0 {
                local lhs = trim(substr(`"`content'"', 1, `=`eqpos'-1'))
                local rhs = trim(substr(`"`content'"', `=`eqpos'+1', .))
                if "`dvar'" == "" & `d_endog' == 0 {
                    local dvar "`lhs'"
                    local ivd  "`rhs'"
                    local d_endog = 1
                }
                else {
                    local mvar "`lhs'"
                    local ivm  "`rhs'"
                }
            }
        }
        else {
            if "`dvar'" == "" & `d_endog' == 0 {
                local dvar "`tok'"
            }
            else {
                local controls "`controls' `tok'"
            }
        }
    }

    *-- validate --
    if "`dvar'" == "" {
        di as err "treatment variable not specified"
        exit 198
    }
    if "`mvar'" == "" | "`ivm'" == "" {
        di as err "mediator equation (m=...) is required"
        exit 198
    }
    if `d_endog' & "`ivd'" == "" {
        di as err "instruments for D missing in (d=...)"
        exit 198
    }
    confirm variable `dvar'
    confirm variable `mvar'

    *========================================================
    * Mark estimation sample
    *========================================================
    marksample touse
    markout `touse' `yvar' `dvar' `mvar' `ivm'
    if `d_endog' markout `touse' `ivd'

    if "`controls'" != "" {
        foreach token of local controls {
            if regexm("`token'", "^i\.") {
                local base = regexr("`token'", "^i\.", "")
            }
            else {
                local base = "`token'"
            }
            markout `touse' `base'
        }
    }

    *-- expand factor variables for coefficient retrieval --
    local ctrl_exp ""
    if "`controls'" != "" {
        qui fvexpand `controls' if `touse'
        local ctrl_exp = r(varlist)
    }

    tempname d_fit m_fit dm_fit Bm bc

    if `d_endog' {
        *========================================================
        * 3-Stage: D endogenous
        *========================================================
        * Stage 1: D on D-instruments + controls -> fitted D
        reg `dvar' `ivd' `controls' if `touse', vce(r)
        predict `d_fit' if `touse', xb

        * Stage 2: M on fitted D + M-instruments + controls
        reg `mvar' `d_fit' `ivm' `controls' if e(sample), vce(r)
        local b0 = _b[_cons]
        local b1 = _b[`d_fit']
        matrix `Bm' = e(b)
        predict `m_fit' if e(sample), xb
        gen `dm_fit' = `d_fit' * `m_fit'

        * Stage 3: Y on fitted D, fitted M, D*M + controls
        reg `yvar' `d_fit' `m_fit' `dm_fit' `controls' if e(sample), vce(r)
        local c1 = _b[`d_fit']
        local c2 = _b[`m_fit']
        local c3 = _b[`dm_fit']
    }
    else {
        *========================================================
        * 2-Stage: D exogenous, M endogenous
        *   D appears in M equation (as in original mediateiv1)
        *   Stage 1: M on D + M-instruments + controls
        *   Stage 2: Y on D, fitted M, D*M + controls
        *========================================================
        reg `mvar' `dvar' `ivm' `controls' if `touse', vce(r)
        local b0 = _b[_cons]
        local b1 = _b[`dvar']
        matrix `Bm' = e(b)
        predict `m_fit' if `touse', xb
        gen `dm_fit' = `dvar' * `m_fit'

        * Stage 2: Y on D, fitted M, D*M + controls
        reg `yvar' `dvar' `m_fit' `dm_fit' `controls' if e(sample), vce(r)
        local c1 = _b[`dvar']
        local c2 = _b[`m_fit']
        local c3 = _b[`dm_fit']
    }

    *-- save N before any r-class commands --
    local nobs = e(N)

    *========================================================
    * S = sum_k ( b_k * mean_k ) over M-equation covariates
    *     (M-instruments + controls incl. factor dummies)
    *========================================================
    local S = 0

    * M-instruments
    foreach v of local ivm {
        scalar `bc' = `Bm'[1, colnumb(`Bm', "`v'")]
        qui su `v' if e(sample)
        local S = `S' + `bc' * r(mean)
    }

    * Controls (factor or continuous)
    if "`ctrl_exp'" != "" {
        qui count if e(sample)
        local ntot = r(N)
        foreach v of local ctrl_exp {
            local col = colnumb(`Bm', "`v'")
            if `col' < . {
                scalar `bc' = `Bm'[1, `col']
            }
            else {
                scalar `bc' = 0
            }
            if regexm("`v'", "^([0-9]+)[boc]*\.") {
                local lvl  = regexs(1)
                local bvar = regexr("`v'", "^([0-9]+)[boc]*\.", "")
                qui count if `bvar' == `lvl' & e(sample)
                local S = `S' + `bc' * (r(N) / `ntot')
            }
            else {
                qui su `v' if e(sample)
                local S = `S' + `bc' * r(mean)
            }
        }
    }

    *========================================================
    * Mediation parameters
    *   NIE(d) = b1 * (c2 + d*c3)
    *   NDE(d) = c1 + (b0 + d*b1)*c3 + c3*S
    *   TE     = NIE(1) + NDE(0)
    *   PM(d)  = NIE(d) / TE
    *========================================================
    local nie0 = `b1' * `c2'
    local nie1 = `b1' * (`c2' + `c3')
    local nde0 = `c1' + `b0'  * `c3' + `c3' * `S'
    local nde1 = `c1' + (`b0' + `b1') * `c3' + `c3' * `S'
    local te   = `nie1' + `nde0'

    if `te' != 0 {
        local pm0 = `nie0' / `te'
        local pm1 = `nie1' / `te'
    }
    else {
        local pm0 = .
        local pm1 = .
    }

    return scalar nie0 = `nie0'
    return scalar nie1 = `nie1'
    return scalar nde0 = `nde0'
    return scalar nde1 = `nde1'
    return scalar te   = `te'
    return scalar pm0  = `pm0'
    return scalar pm1  = `pm1'

    *========================================================
    * Output
    *========================================================
    local mtitle "3SLS"
    if !`d_endog' local mtitle "2SLS"

    if `reps' > 0 {
        *-- preserve if/in for recursive bootstrap call --
        local ifin ""
        if `"`if'"' != "" local ifin `"`if'"'
        if `"`in'"' != "" local ifin `"`ifin' `in'"'

        if "`seed'" != "" set seed `seed'

        *-- Save key values as scalars before bootstrap --
        *   (local macros may not survive the bootstrap prefix)
        scalar _medv_dendog = `d_endog'
        scalar _medv_nobs   = `nobs'
        scalar _medv_nreps  = `reps'

        if `d_endog' {
            bootstrap nie0=r(nie0) nie1=r(nie1) nde0=r(nde0) nde1=r(nde1) ///
                te=r(te) pm0=r(pm0) pm1=r(pm1), reps(`reps') noheader notable: ///
                mediateiv `yvar' (`dvar'=`ivd') (`mvar'=`ivm') `controls' `ifin', reps(0)
        }
        else {
            bootstrap nie0=r(nie0) nie1=r(nie1) nde0=r(nde0) nde1=r(nde1) ///
                te=r(te) pm0=r(pm0) pm1=r(pm1), reps(`reps') noheader notable: ///
                mediateiv `yvar' `dvar' (`mvar'=`ivm') `controls' `ifin', reps(0)
        }

        *-- Retrieve saved values from scalars --
        local _dendog = scalar(_medv_dendog)
        local _nobs   = scalar(_medv_nobs)
        local _nreps  = scalar(_medv_nreps)
        scalar drop _medv_dendog _medv_nobs _medv_nreps

        *-- Bootstrap table (inline, no subroutine) --
        tempname Bs Vs
        matrix `Bs' = e(b)
        matrix `Vs' = e(V)

        local effects "nie0 nie1 nde0 nde1 te pm0 pm1"
        local labels `""NIE (D=0)" "NIE (D=1)" "NDE (D=0)" "NDE (D=1)" "TE" "PM (D=0)" "PM (D=1)""'

        local _mtitle "3SLS"
        if !`_dendog' local _mtitle "2SLS"

        di _n as txt "IV Mediation Analysis (`_mtitle' with D*M Interaction)"
        di as txt "Bootstrap replications = " as res `_nreps'
        di as txt "{hline 72}"
        di as txt _col(2)  "Effect" ///
           _col(16) "Estimate" ///
           _col(30) "Std. Err." ///
           _col(44) "[95% Conf. Interval]" ///
           _col(68) "P>|z|"
        di as txt "{hline 72}"

        local i = 1
        foreach ef of local effects {
            local lbl : word `i' of `labels'
            local b   = `Bs'[1, `i']
            local se  = sqrt(`Vs'[`i', `i'])
            local z   = `b' / `se'
            local p   = 2 * (1 - normal(abs(`z')))
            local lo  = `b' - 1.96 * `se'
            local hi  = `b' + 1.96 * `se'

            di as txt _col(2) "`lbl'" ///
               as res _col(16) %10.7f `b'  ///
               _col(30) %10.7f `se' ///
               _col(44) %10.7f `lo' ///
               _col(56) %10.7f `hi' ///
               _col(70) %5.3f `p'
            local i = `i' + 1
        }
        di as txt "{hline 72}"
        di as txt "Observations" _col(16) as res %10.0fc `_nobs'
    }
    else {
        di _n as txt "IV Mediation Analysis (`mtitle' with D*M Interaction)"
        di as txt "{hline 60}"
        di as txt _col(6) "Effect" _col(28) "Estimate"
        di as txt "{hline 60}"
        di as txt _col(6) "NIE (D=0)"  _col(28) as res %10.7f `nie0'
        di as txt _col(6) "NIE (D=1)"  _col(28) as res %10.7f `nie1'
        di as txt _col(6) "NDE (D=0)"  _col(28) as res %10.7f `nde0'
        di as txt _col(6) "NDE (D=1)"  _col(28) as res %10.7f `nde1'
        di as txt _col(6) "TE"         _col(28) as res %10.7f `te'
        di as txt _col(6) "PM (D=0)"   _col(28) as res %10.7f `pm0'
        di as txt _col(6) "PM (D=1)"   _col(28) as res %10.7f `pm1'
        di as txt "{hline 60}"
        di as txt "Observations" _col(28) as res %10.0fc `nobs'
    }
end