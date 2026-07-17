{smcl}
{* *! mediateiv.sthlp 16jul2026}{...}
{title:mediateiv}

{p2colset 5 18 22 2}{...}
{p2col:{cmd:mediateiv} {hline 2}}IV causal mediation analysis with treatment-mediator interaction{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{pstd}
When the treatment {it:D} is endogenous (instrumented):

{p 8 15 2}
{cmd:mediateiv} {depvar} {cmd:(}{it:dvar} {cmd:=} {it:ivd_list}{cmd:)} {cmd:(}{it:mvar} {cmd:=} {it:ivm_list}{cmd:)} [{it:controls}] {ifin}
{cmd:,} [{opth reps:(integer)} {opth seed:(string)}]

{pstd}
When the treatment {it:D} is exogenous:

{p 8 15 2}
{cmd:mediateiv} {depvar} {it:dvar} {cmd:(}{it:mvar} {cmd:=} {it:ivm_list}{cmd:)} [{it:controls}] {ifin}
{cmd:,} [{opth reps:(integer)} {opth seed:(string)}]

{marker syntax_element}{...}
{synoptset 22 tabbed}{...}
{synopthdr:element}
{synoptline}
{synopt:{depvar}}dependent variable (outcome){p_end}
{synopt:{cmd:(}{it:dvar} {cmd:=} {it:ivd_list}{cmd:)}}treatment equation: {it:dvar} is the treatment variable, {it:ivd_list} are instruments for {it:dvar}; presence of this equation triggers 3-stage estimation{p_end}
{synopt:{it:dvar}}treatment variable specified directly (no instruments); triggers 2-stage estimation{p_end}
{synopt:{cmd:(}{it:mvar} {cmd:=} {it:ivm_list}{cmd:)}}mediator equation: {it:mvar} is the mediator, {it:ivm_list} are instruments for {it:mvar}{p_end}
{synopt:{it:controls}}covariates (controls); {cmd:i.} factor-variable notation is supported, e.g., {cmd:i.prov}{p_end}
{synoptline}
{p2colreset}{...}

{marker options}{...}
{title:Options}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opth reps:(integer)}}Number of bootstrap replications. Default is {cmd:reps(100)}. {cmd:reps(0)} suppresses bootstrap and reports point estimates only. {cmd:reps(500)} or higher is recommended for stable inference.{p_end}
{synopt:{opth seed:(string)}}Random-number seed for reproducible bootstrap results, e.g., {cmd:seed(12345)}.{p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
{cmd:mediateiv} performs instrumental-variable (IV) causal mediation analysis with a treatment-mediator interaction ({it:D}*{it:M}). It decomposes the total effect (TE) of the treatment {it:D} on the outcome {it:Y} into a natural indirect effect (NIE, through the mediator {it:M}) and a natural direct effect (NDE), while addressing endogeneity of the treatment and/or the mediator via instruments.

{pstd}
Two estimation modes are supported:

{pstd}
{bf:1. Three-stage estimation (D endogenous).} Use the syntax {cmd:(}{it:dvar}{cmd:=}{it:ivd_list}{cmd:)} when the treatment is endogenous.

{p 8 12 2}
{it:Stage 1:} {cmd:reg} {it:dvar ivd_list controls}{cmd:, vce(r)} -> fitted D

{p 8 12 2}
{it:Stage 2:} {cmd:reg} {it:mvar fitted_D ivm_list controls}{cmd:, vce(r)} -> fitted M, coefficients b0, b1

{p 8 12 2}
{it:Stage 3:} {cmd:reg} {it:yvar fitted_D fitted_M fitted_D*fitted_M controls}{cmd:, vce(r)} -> coefficients c1, c2, c3

{pstd}
{bf:2. Two-stage estimation (D exogenous).} Use the syntax {it:dvar} (without parentheses) when the treatment is exogenous. The treatment still appears in the mediator equation.

{p 8 12 2}
{it:Stage 1:} {cmd:reg} {it:mvar dvar ivm_list controls}{cmd:, vce(r)} -> fitted M, coefficients b0, b1

{p 8 12 2}
{it:Stage 2:} {cmd:reg} {it:yvar dvar fitted_M dvar*fitted_M controls}{cmd:, vce(r)} -> coefficients c1, c2, c3

{pstd}
All regressions use heteroskedasticity-robust standard errors ({cmd:vce(r)}).

{marker formulas}{...}
{title:Mediation parameters}

{pstd}
Let {it:b0} = intercept, {it:b1} = coefficient on D (or fitted D) in the M equation; {it:c1} = coefficient on D, {it:c2} = coefficient on M, {it:c3} = coefficient on {it:D}*{it:M} in the Y equation; and {it:S} = {it:sum_k(b_k * mean_k)} over the covariates in the M equation (for factor variables, the group proportion is used as the mean). Then:

{p 8 12 2}
{it:NIE(d)} = {it:b1} * ({it:c2} + {it:d} * {it:c3})

{p 8 12 2}
{it:NDE(d)} = {it:c1} + ({it:b0} + {it:d} * {it:b1}) * {it:c3} + {it:c3} * {it:S}

{p 8 12 2}
{it:TE} = {it:NIE(1)} + {it:NDE(0)}

{p 8 12 2}
{it:PM(d)} = {it:NIE(d)} / {it:TE}

{pstd}
where {it:d} = 0 or 1 indexes the treatment level at which the effect is evaluated.

{marker output}{...}
{title:Output}

{pstd}
{bf:Point estimates} ({cmd:reps(0)}): a table reporting NIE(D=0), NIE(D=1), NDE(D=0), NDE(D=1), TE, PM(D=0), and PM(D=1), with sample size at the bottom.

{pstd}
{bf:Bootstrap inference} ({cmd:reps(}#{cmd:)} > 0, the default being 100): a table reporting each effect's estimate, bootstrap standard error, 95% confidence interval, and {it:P>|z|}, with sample size at the bottom.

{marker saved_results}{...}
{title:Saved results}

{pstd}
{cmd:mediateiv} stores the following in {cmd:r()}:

{synoptset 15 tabbed}{...}
{synopthdr:scalars}
{synoptline}
{synopt:{cmd:r(nie0)}}natural indirect effect at D=0{p_end}
{synopt:{cmd:r(nie1)}}natural indirect effect at D=1{p_end}
{synopt:{cmd:r(nde0)}}natural direct effect at D=0{p_end}
{synopt:{cmd:r(nde1)}}natural direct effect at D=1{p_end}
{synopt:{cmd:r(te)}}total effect{p_end}
{synopt:{cmd:r(pm0)}}proportion mediated at D=0{p_end}
{synopt:{cmd:r(pm1)}}proportion mediated at D=1{p_end}
{synoptline}
{p2colreset}{...}

{pstd}
When {cmd:reps(}#{cmd:)} > 0, {cmd:e(b)} and {cmd:e(V)} from {helpb bootstrap} are also available, containing the bootstrap estimates and variance-covariance matrix of the seven mediation parameters (in the order {it:nie0, nie1, nde0, nde1, te, pm0, pm1}).

{marker examples}{...}
{title:Examples}

{pstd}
{bf:Setup}

{phang2}{cmd:. use urban0, clear}

{phang2}{cmd:. global xlist age age2 i.prov}

{pstd}
{bf:D endogenous — three-stage estimation with bootstrap}

{phang2}{cmd:. mediateiv lnwage (college=focup moccup) (occup=focup moccup) $xlist, reps(500) seed(12345)}

{pstd}
{bf:D exogenous — two-stage estimation with bootstrap}

{phang2}{cmd:. mediateiv lnwage college (occup=focup moccup) $xlist, reps(500) seed(12345)}

{pstd}
{bf:Point estimates only (no bootstrap)}

{phang2}{cmd:. mediateiv lnwage (college=focup moccup) (occup=focup moccup) $xlist}

{phang2}{cmd:. mediateiv lnwage college (occup=focup moccup) $xlist}

{marker remarks}{...}
{title:Remarks}

{pstd}
1. The treatment-mediator interaction {it:D}*{it:M} is always included in the outcome equation. If {it:c3} = 0 (no interaction), NIE(0) = NIE(1) and NDE(0) = NDE(1).

{pstd}
2. Factor-variable notation (e.g., {cmd:i.prov}) is supported in {it:controls}. The underlying categories are automatically expanded for coefficient retrieval, and group proportions are used as means in the {it:S} calculation.

{pstd}
3. When {cmd:reps(}#{cmd:)} > 0, bootstrap is implemented via recursive calls to {cmd:mediateiv} with {cmd:reps(0)} to obtain point estimates for each bootstrap sample. This ensures the same estimation code is used for both bootstrap replicates and point estimates.

{pstd}
4. Bootstrap inference uses a normal approximation: {it:z} = estimate / SE, {it:p} = 2*(1 - Phi(|z|)), and 95% CI = estimate +/- 1.96*SE.

{marker author}{...}
{title:Author}

{pstd}
Xiliang Zhao, Xiamen University{break}
zhaoxiliang@gmail.com

{marker also}{...}
{title:Also see}

{psee}
Help: {helpb regress}, {helpb bootstrap}, {helpb ivregress}
{p_end}
