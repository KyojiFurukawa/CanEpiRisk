# Calculate Years of Life Lost (YLL) due to Radiation Exposure

Computes the expected **Years of Life Lost (YLL)** attributable to
radiation exposure. This function integrates site-specific baseline
mortality and all-cause mortality rates with user-specified excess risk
models (ERR/EAR) to estimate the expected reduction in life expectancy
due to radiation-attributable cancer mortality.

The function supports both single-age and protracted exposure scenarios
and allows for uncertainty quantification via Monte Carlo sampling.

## Usage

``` r
YLL(exposure, reference, riskmodel, option)
```

## Arguments

- exposure:

  A **list** specifying the exposure scenario:

  - `agex` — age(s) at exposure (scalar or numeric vector).

  - `doseGy` — dose(s) in Gray (Gy), scalar or vector of same length as
    `agex`.

  - `sex` — sex indicator (1 = male, 2 = female).

- reference:

  A **list** specifying baseline and all-cause mortality information:

  - `baseline` — data.frame of site-specific mortality (or incidence)
    rates per person-year, including columns `age`, `male`, and
    `female`.

  - `mortality` — data.frame of all-cause mortality rates (per
    person-year) with same structure and age range as `baseline`.

  These should correspond to the same population or region.

- riskmodel:

  A **list** defining the excess risk model, typically of the form:

  - `err` and/or `ear` — sublists each containing:

    - `para` — vector of model parameters.

    - `var` — variance–covariance matrix of parameters (for
      multi-parameter models).

    - `ci` — 95% confidence bounds for a single-parameter model
      (alternative to `var`).

    - `f(beta, data, lag)` — function returning age-specific excess
      risks.

- option:

  A **list** of optional settings:

  - `maxage` — maximum attained age for accumulation (default: 100).

  - `err_wgt` — weight between ERR (1) and EAR (0) transfer models.

  - `n_mcsamp` — number of Monte Carlo draws for uncertainty estimation
    (default: 10,000).

  - `alpha` — significance level for confidence intervals (default:
    0.05).

## Value

A named numeric vector with point and interval summaries of cumulative
excess risk, typically including:

- `mle`: point estimate,

- `mean`, `median`: Monte Carlo summaries,

- `ci_lo`, `ci_up`: confidence interval.

Values are per person; multiply by `1e4` or `1e5` to report per 10,000
or 100,000.

## Details

The YLL estimate integrates the excess mortality over age, weighted by
remaining life expectancy, producing an interpretable measure of life
expectancy loss attributable to radiation. The computation uses region-
and sex-specific mortality data provided in the package or supplied by
the user.

Typical usage involves pairing `YLL()` with
[`CER()`](https://kyojifurukawa.github.io/CanEpiRisk/reference/CER.md)
or
[`population_YLL()`](https://kyojifurukawa.github.io/CanEpiRisk/reference/population_YLL.md)
for complementary assessments of probability-based and life loss–based
risk metrics.

## See also

[`CER`](https://kyojifurukawa.github.io/CanEpiRisk/reference/CER.md),
[`population_YLL`](https://kyojifurukawa.github.io/CanEpiRisk/reference/population_YLL.md),
[`population_LAR`](https://kyojifurukawa.github.io/CanEpiRisk/reference/population_LAR.md)

## Examples

``` r
set.seed(100)
## Example 1: All-solid cancer mortality (Region 1, female, 0.1 Gy at age 15)
exp1 <- list(agex = 15, doseGy = 0.1, sex = 2)
ref1 <- list(
  baseline  = Mortality[[1]]$allsolid,
  mortality = Mortality[[1]]$allcause
)
mod1 <- LSS_mortality$allsolid$L
opt1 <- list(maxage = 100, err_wgt = 1, n_mcsamp = 10000)

YLL(exposure = exp1, reference = ref1, riskmodel = mod1, option = opt1)
#>         mle        mean      median  ci_lo.2.5% ci_up.97.5% 
#>   0.2408580   0.2437128   0.2417327   0.1884222   0.3090655 

## Example 2: Leukaemia incidence (Region 4), male, 100 mGy evenly across ages 30–45
exp2 <- list(agex = 30:44 + 0.5, doseGy = rep(0.1 / 15, 15), sex = 1)
ref2 <- list(
  baseline  = Incidence[[4]]$leukaemia,
  mortality = Mortality[[4]]$allcause
)
mod2 <- LSS_incidence$leukaemia$LQ
opt2 <- list(maxage = 60, err_wgt = 0, n_mcsamp = 10000)

YLL(exposure = exp2, reference = ref2, riskmodel = mod2, option = opt2)
#>           mle          mean        median    ci_lo.2.5%   ci_up.97.5% 
#>  4.118829e-03  4.046922e-03  4.046761e-03 -6.836959e-05  8.171657e-03 
```
