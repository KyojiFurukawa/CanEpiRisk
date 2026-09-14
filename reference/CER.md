# CER: Calculating Cumulative Excess Risks

Compute the **cumulative excess risk (CER)** attributable to radiation
exposure for an individual scenario or cohort, integrating excess risks
over age up to a specified maximum while accounting for competing
all-cause mortality and model uncertainty. This is the core engine for
lifetime risk estimation in **CanEpiRisk**.

## Usage

``` r
CER(exposure, reference, riskmodel, option)
```

## Arguments

- exposure:

  list. Exposure scenario with components:

  - `agex`: age(s) at exposure (scalar or vector).

  - `doseGy`: dose(s) in gray (Gy); same length as `agex` if vectorized.

  - `sex`: sex indicator (`1` = male, `2` = female).

- reference:

  list. Baseline reference data for the same population/region with:

  - `baseline`: data frame of site-specific baseline rates (incidence or
    mortality), with columns `age`, `male`, `female` on ages 1:100.

  - `mortality`: data frame of all-cause mortality (same columns/age
    grid).

- riskmodel:

  list. Radiation risk model definition with sublists for *excess
  relative risk* (ERR) and *excess absolute risk* (EAR), e.g.:

  - `err`/`ear`: each a list containing `para` (numeric parameter
    vector), `var` (variance–covariance matrix) *or* `ci` (confidence
    interval for 1-parameter models), and `f` (function of the form
    `f(beta, data, lag)` returning age-specific excess risk).

- option:

  list. Optional settings:

  - `maxage`: maximum attained age for accumulation (e.g., `100`).

  - `err_wgt`: weight to blend ERR vs EAR (`1` = pure ERR; `0` = pure
    EAR; intermediate values allowed).

  - `n_mcsamp`: Monte Carlo sample size for uncertainty propagation
    (e.g., `10000`).

  - `alpha`: significance level for interval estimation (default
    `0.05`).

## Value

A named numeric vector with point and interval summaries of cumulative
excess risk, typically including:

- `mle`: point estimate,

- `mean`, `median`: Monte Carlo summaries,

- `ci_lo`, `ci_up`: confidence intervals.

Values are per person; multiply by `1e4` or `1e5` to report per 10,000
or 100,000.

## Details

Calculating Cumulative Excess Risks (CER)

The function integrates age-specific excess risks implied by the
supplied ERR/EAR model(s) across the attained-age grid (typically ages
1–100). When `n_mcsamp > 1`, parameter uncertainty is propagated via
Monte Carlo sampling using either the model variance–covariance matrix
(`var`) or a 95% CI (`ci`) for 1-parameter models. For protracted
exposures (vectorized `agex`/`doseGy`), contributions are summed across
exposure segments. Baseline site-specific rates and all-cause mortality
must correspond to the same population/region and share the same age
grid and sex coding.

## Units & Alignment

Doses must be in Gy. Ensure `baseline` and `mortality` tables are from
the same population/region and aligned on age (1–100) and sex coding.

## See also

[`population_LAR`](https://kyojifurukawa.github.io/CanEpiRisk/reference/population_LAR.md),
[`population_YLL`](https://kyojifurukawa.github.io/CanEpiRisk/reference/population_YLL.md),
[`YLL`](https://kyojifurukawa.github.io/CanEpiRisk/reference/YLL.md)

## Examples

``` r
set.seed(100)
# Example 1: All-solid cancer mortality (Region 1), female, 0.1 Gy at age 15
exp1 <- list(agex = 15, doseGy = 0.1, sex = 2)
ref1 <- list(
  baseline  = Mortality[[1]]$allsolid,
  mortality = Mortality[[1]]$allcause
)
mod1 <- LSS_mortality$allsolid$L
opt1 <- list(maxage = 100, err_wgt = 1, n_mcsamp = 10000)

CER(exposure = exp1, reference = ref1, riskmodel = mod1, option = opt1) * 10000
#>         mle        mean      median  ci_lo.2.5% ci_up.97.5% 
#>    156.0667    158.3529    156.4740    114.7944    214.0962 

# Example 2: Leukaemia incidence (Region 4), male, 100 mGy evenly across ages 30–45
exp2 <- list(agex = 30:44 + 0.5, doseGy = rep(0.1/15, 15), sex = 1)
ref2 <- list(
  baseline  = Incidence[[4]]$leukaemia,
  mortality = Mortality[[4]]$allcause
)
mod2 <- LSS_incidence$leukaemia$LQ
opt2 <- list(maxage = 60, err_wgt = 0, n_mcsamp = 10000)

CER(exposure = exp2, reference = ref2, riskmodel = mod2, option = opt2) * 10000
#>         mle        mean      median  ci_lo.2.5% ci_up.97.5% 
#>  3.68374227  3.62558094  3.63197695 -0.06062042  7.27448084 

```
