# Generating s Monte Carlo sample of CER

Generatie s Monte Carlo sample of CER from a risk model under a
specified exposure scenario.

## Usage

``` r
mc_CER(exposure, reference, riskmodel, option)
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

information of calculated excess risk (data.frame)

## Examples

``` r
 # The following examples use default data provided in CanEpiRisk package
 # for riskmodels (LSS_mortality and LSS_incidence) derived from Life Span Study
 # and baseline mortality and incidence rates for WHO global regions (Mortality and Incidence).

 # Example 1: allsolid mortality, Region-1, female, 0.1Gy at age 15, followed up to age 100, LSS linear ERR
 exp1 <- list( agex=5, doseGy=0.1, sex=2 )   # exposure scenario
 ref1 <- list( baseline=Mortality[[1]]$allsolid,        # baseline rates
              mortality=Mortality[[1]]$allcause )       # all-cause mortality
 mod1 <- LSS_mortality$allsolid$L                       # risk model
 opt1 <- list( maxage=100, err_wgt=1, n_mcsamp=10000 )  # option
 CER(  exposure=exp1, reference=ref1, riskmodel=mod1, option=opt1 ) * 10000 # cases per 10,000
#>         mle        mean      median  ci_lo.2.5% ci_up.97.5% 
#>    221.1383    226.9980    221.3491    147.4071    337.7613 

```
