# mc_CER: Calculating excess risks

Calculate the excess risk from a risk model under a specified exposure
scenario.

## Usage

``` r
Comp_Exrisk(exposure, riskmodel, option, per = 1)
```

## Arguments

- exposure:

  a list object that specifies the exposure scenario, which contains
  `agex` (a single value or a vector for age(s) at exposure), 'doseGy'
  (a single value or a vector of dose(s) in Gy), and 'sex' (1 or 2 for
  male or female).

- riskmodel:

  a list object that specifies the risk model, which contains two list
  objects named `err` for excess relative rate model and 'ear' for
  excess absolute rate model, each of which contains a vector 'para' for
  model parameter estimates and a function 'f' to compute the excess
  risk given a parameter vector and exposure information (e.g., dose,
  age at exposure, sex, attained age).

- option:

  a list object that specifies optional settings for risk calculation,
  which contains an integer value 'maxage' for the maximum age to follow
  up and a value 'err_wgt' for the weight for risk transfer (1=err,
  0=ear).

- per:

  an integer value for the risk denominator (default=1).

## Value

information of calculated excess risk (data.frame)

## See also

[LSS_mortality](https://kyojifurukawa.github.io/CanEpiRisk/reference/LSS_mortality.md),
[LSS_incidence](https://kyojifurukawa.github.io/CanEpiRisk/reference/LSS_incidence.md)

## Examples

``` r
 # The following examples use default data provided in CanEpiRisk package
 # for riskmodels (LSS_mortality and LSS_incidence) derived from Life Span Study
 # and baseline mortality and incidence rates for WHO global regions (Mortality and Incidence).

```
