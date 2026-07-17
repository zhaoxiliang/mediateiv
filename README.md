# mediateiv

**IV Causal Mediation Analysis with Treatment–Mediator Interaction for Stata**

## Overview

`mediateiv` performs instrumental-variable (IV) causal mediation analysis with a treatment–mediator interaction (D×M), which developed by Zhao(2025). It decomposes the total effect (TE) of a treatment *D* on an outcome *Y* into a natural indirect effect (NIE, operating through the mediator *M*) and a natural direct effect (NDE), while addressing endogeneity of the treatment and/or the mediator via instruments.

Two estimation modes are supported:

| Mode | When to use | Syntax |
|------|-------------|--------|
| **3-stage (3SLS)** | Treatment *D* is endogenous | `mediateiv y (d=z1) (m=z2) controls` |
| **2-stage (2SLS)** | Treatment *D* is exogenous | `mediateiv y d (m=z2) controls` |

All regressions use heteroskedasticity-robust standard errors (`vce(r)`). Bootstrap inference is built in for standard errors, confidence intervals, and *p*-values.

## Installation

Copy `mediateiv.ado` and `mediateiv.sthlp` to the directory of Stata "\ado\plus" or "\ado\personal", then you can use it as native Stata command.

Then verify:

```stata
which mediateiv
help mediateiv
```

## Syntax

### D endogenous — three-stage estimation

```
mediateiv depvar (dvar = ivd_list) (mvar = ivm_list) [controls] [if] [in] [, reps(#) seed(#)]
```

### D exogenous — two-stage estimation

```
mediateiv depvar dvar (mvar = ivm_list) [controls] [if] [in] [, reps(#) seed(#)]
```

### Elements

| Element | Description |
|---------|-------------|
| `depvar` | Dependent variable (outcome) |
| `(dvar = ivd_list)` | Treatment equation: `dvar` is the treatment; `ivd_list` are instruments for `dvar`. Presence of this equation triggers 3-stage estimation. |
| `dvar` | Treatment variable specified directly (no instruments). Triggers 2-stage estimation. |
| `(mvar = ivm_list)` | Mediator equation: `mvar` is the mediator; `ivm_list` are instruments for `mvar`. |
| `controls` | Covariates. `i.` factor-variable notation is supported (e.g., `i.prov`). |
| `if` / `in` | Standard Stata qualifiers for restricting the estimation sample. |

### Options

| Option | Description |
|--------|-------------|
| `reps(#)` | Number of bootstrap replications. Default is `reps(100)`. `reps(0)` suppresses bootstrap and reports point estimates only. `reps(500)` or higher is recommended for stable inference. |
| `seed(#)` | Random-number seed for reproducible bootstrap results. |

## Method

### Three-stage estimation (D endogenous)

| Stage | Regression | Output |
|-------|-----------|--------|
| 1 | `reg dvar ivd_list controls, vce(r)` | Fitted D |
| 2 | `reg mvar fitted_D ivm_list controls, vce(r)` | Fitted M; coefficients b₀, b₁ |
| 3 | `reg yvar fitted_D fitted_M fitted_D×fitted_M controls, vce(r)` | Coefficients c₁, c₂, c₃ |

### Two-stage estimation (D exogenous)

| Stage | Regression | Output |
|-------|-----------|--------|
| 1 | `reg mvar dvar ivm_list controls, vce(r)` | Fitted M; coefficients b₀, b₁ |
| 2 | `reg yvar dvar fitted_M dvar×fitted_M controls, vce(r)` | Coefficients c₁, c₂, c₃ |

The only difference between the two modes is whether D is instrumented (3-stage) or used directly (2-stage). In both cases D appears in the mediator equation.

### Mediation parameters

Let **b₀** = intercept, **b₁** = coefficient on D in the M equation; **c₁** = coefficient on D, **c₂** = coefficient on M, **c₃** = coefficient on D×M in the Y equation; and **S** = Σ(b_k × mean_k) over the covariates in the M equation (for factor variables, the group proportion is used as the mean). Then:

| Parameter | Formula |
|-----------|---------|
| NIE(d) | b₁ × (c₂ + d × c₃) |
| NDE(d) | c₁ + (b₀ + d × b₁) × c₃ + c₃ × S |
| TE | NIE(1) + NDE(0) |
| PM(d) | NIE(d) / TE |

where *d* = 0 or 1 indexes the treatment level at which the effect is evaluated.

## Saved results

`mediateiv` stores the following in `r()`:

| Scalar | Description |
|--------|-------------|
| `r(nie0)` | Natural indirect effect at D = 0 |
| `r(nie1)` | Natural indirect effect at D = 1 |
| `r(nde0)` | Natural direct effect at D = 0 |
| `r(nde1)` | Natural direct effect at D = 1 |
| `r(te)` | Total effect |
| `r(pm0)` | Proportion mediated at D = 0 |
| `r(pm1)` | Proportion mediated at D = 1 |

When `reps(#) > 0`, `e(b)` and `e(V)` from `bootstrap` are also available, containing the bootstrap estimates and variance–covariance matrix of the seven mediation parameters (in order: nie0, nie1, nde0, nde1, te, pm0, pm1).

## Examples

```stata
* Load data and prepare sample
use urban0, clear
keep if edu > 3
qui reg linc college occup seduc age age2 male married i.prov moccup foccup
keep if e(sample)

* --- D endogenous: three-stage estimation with bootstrap ---
mediateiv linc (college=seduc) (occup=focup moccup) age age2 male married i.prov, reps(500) seed(12345)

* --- D exogenous: two-stage estimation with bootstrap ---
mediateiv linc college (occup=focup moccup) age age2 male married i.prov, reps(500) seed(12345)

* --- Point estimates only (no bootstrap) ---
mediateiv linc (college=seduc) (occup=focup moccup) age age2 male married i.prov, reps(0)
```

## Sample output

**Bootstrap inference** (`reps(500)`):

```
IV Mediation Analysis (3SLS with D*M Interaction)
Bootstrap replications = 500
------------------------------------------------------------------------
  Effect         Estimate    Std. Err.   [95% Conf. Interval]   P>|z|
------------------------------------------------------------------------
  NIE (D=0)     0.0739770   0.0672757    -0.0579487  0.2059027  0.272
  NIE (D=1)     0.2570077   0.0340070     0.1903447  0.3236708  0.000
  NDE (D=0)     0.3810872   0.0339541     0.3145337  0.4476407  0.000
  NDE (D=1)     0.5641178   0.0574278     0.4515360  0.6766995  0.000
  TE            0.6380948   0.0349941     0.5695107  0.7066790  0.000
  PM (D=0)      0.1159342   0.1046832    -0.0892721  0.3211405  0.268
  PM (D=1)      0.4027735   0.0497866     0.3051893  0.5003577  0.000
------------------------------------------------------------------------
Observations        3139
```

**Point estimates only** (`reps(0)`):

```
IV Mediation Analysis (3SLS with D*M Interaction)
------------------------------------------------------------
     Effect                  Estimate
------------------------------------------------------------
     NIE (D=0)              0.0739770
     NIE (D=1)              0.2570077
     NDE (D=0)              0.3810872
     NDE (D=1)              0.5641178
     TE                     0.6380948
     PM (D=0)               0.1159342
     PM (D=1)               0.4027735
------------------------------------------------------------
Observations                  3139
```

## Files

| File | Description |
|------|-------------|
| `mediateiv.ado` | Main program file |
| `mediateiv.sthlp` | Stata help file (`help mediateiv`) |
| `urban0.dta` | Example dataset (CHIP2018 urban households) |

## Remarks

1. The treatment–mediator interaction D×M is always included in the outcome equation. If c₃ = 0 (no interaction), NIE(0) = NIE(1) and NDE(0) = NDE(1).

2. Factor-variable notation (e.g., `i.prov`) is supported in `controls`. Categories are automatically expanded for coefficient retrieval, and group proportions serve as means in the *S* calculation.

3. Bootstrap is implemented via recursive calls to `mediateiv` with `reps(0)` for each resample, ensuring identical estimation code across replicates and point estimates.

4. Bootstrap inference uses a normal approximation: *z* = estimate / SE, *p* = 2 × (1 − Φ(|z|)), 95% CI = estimate ± 1.96 × SE.

## Author

**Xiliang Zhao**, Xiamen University  
Email: zhaoxiliang@gmail.com

## Reference
赵西亮(2025). 因果中介分析的理论进展及其应用. 数量经济技术经济研究，第2期。
## Also see

Stata help: [`regress`](https://www.stata.com/help.cgi?regress), [`bootstrap`](https://www.stata.com/help.cgi?bootstrap), [`ivregress`](https://www.stata.com/help.cgi?ivregress)

## License

This project is released under the MIT License.
