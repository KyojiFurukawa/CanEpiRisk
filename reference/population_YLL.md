# Calculating population-averaged years of life lost

Compute **population-averaged years of life lost (YLL)** due to
radiation exposure, aggregating over a population age distribution and
sex-specific baseline rates. The function combines user-specified excess
risk models (ERR/EAR) with site-specific baseline incidence/mortality
and all-cause mortality to produce age-category summaries and overall
totals, with uncertainty from Monte Carlo sampling.

## Usage

``` r
population_YLL(dsGy, reference, riskmodel, agex = 1:8 * 10 - 5, nmc = 10000)
```

## Arguments

- dsGy:

  Numeric scalar. Radiation **dose** in Gy (or Sv). Must be nonnegative.

- reference:

  List describing the target **reference population**, with components:

  - `baseline` — data frame of site-specific baseline **incidence or
    mortality** rates (per person-year) over ages 1–100, with columns
    `age`, `male`, `female`.

  - `mortality` — data frame of **all-cause mortality** (same
    structure/age grid).

  - `agedist` — data frame or vector with the **population age
    distribution** used to average risks across ages (e.g., by single
    year or grouped ages).

  All three components must refer to the **same region/population** and
  share consistent sex coding.

- riskmodel:

  List defining the radiation **risk model**. Typically contains
  sublists for ERR and/or EAR:

  - `err` / `ear` — each is a list with:

    - `para` — numeric vector of model parameters;

    - `var` — parameter variance–covariance matrix (for multi-parameter
      models), *or* `ci` — confidence interval for one-parameter models;

    - `f` — function of the form `f(beta, data, lag)` returning
      age-specific excess risk given parameters and data (dose, ages,
      sex).

  See package examples (e.g., `LSS_mortality$allsolid$L`,
  `LSS_incidence$leukaemia$LQ`).

- agex:

  Numeric vector of **ages at exposure** (midpoints for grouped
  categories). Default is `1:8 * 10 - 5` (i.e., 5, 15, …, 75)
  representing 0–10, 10–20, …, 70–80.

- nmc:

  Integer Monte Carlo sample size for uncertainty propagation. Default
  `10000`.

## Value

A **list** with components for each transfer model present (typically
`$err` and `$ear`). Each component is a data frame with one row per
exposure age category (rows named by `agex`) plus an `all` row, and
columns:

- `male`, `male_lo`, `male_up`

- `female`, `female_lo`, `female_up`

- `all`, `all_lo`, `all_up`

Values are YLL.

## Details

Let `agedist` denote the population age distribution (grouped ages) and
`baseline`/`mortality` denote site-specific baseline rates and all-cause
mortality for the same region/population and sex coding. For each
age-at-exposure category in `agex`, `population_YLL()` evaluates YLL
under the supplied risk model(s) and then averages across the age
structure, returning sex-specific estimates and an `all`-sex average.
Uncertainty is obtained by drawing model parameters either from a
variance–covariance matrix (`var`) or, in one-parameter models, from 95%
confidence bounds (`ci`) when provided.

The default `agex = 1:8 * 10 - 5` corresponds to midpoint ages 5, 15, …,
75, i.e., exposure categories 0–10, 10–20, …, 70–80 years, which should
match the categories of `agedist`.

## Units & Alignment

Doses must be in Gy (or consistent with your model). Ensure `baseline`
and `mortality` share the same age grid (typically 1–100) and sex
coding, and that `agedist` corresponds to the same population.

## See also

[`YLL`](https://kyojifurukawa.github.io/CanEpiRisk/reference/YLL.md),
[`plot_agedist`](https://kyojifurukawa.github.io/CanEpiRisk/reference/plot_agedist.md),
`population_CER`,
[`CER`](https://kyojifurukawa.github.io/CanEpiRisk/reference/CER.md)

## Examples

``` r
set.seed(100)
# The following examples use default data provided in CanEpiRisk package
# for riskmodels (LSS_mortality) derived from Life Span Study
#     baseline rates and age distribution for WHO riskmodels (Mortality, agedist_rgn)
# Example: allsolid mortality, Region-1, exposed to 0.1 Gy, followed up to age 100, LSS linear ERR
ref1 <- list(  baseline=Mortality[[1]]$allsolid,     # baseline rates
mortality=Mortality[[1]]$allcause,     # allcause mortality
agedist=agedist_rgn[[1]] )           # age distribution
mod1 <- LSS_mortality$allsolid$L                     # risk model
population_YLL( dsGy=0.1, reference=ref1, riskmodel=mod1 )
#> $err
#>            male     male_lo     male_up      female   female_lo  female_up
#> 5   0.172093332 0.115054187 0.251928543 0.347325544 0.246518841 0.48922054
#> 15  0.120103006 0.084205848 0.164594216 0.241732716 0.188422239 0.30906550
#> 25  0.083974983 0.058456033 0.115108413 0.168027734 0.134188358 0.20943433
#> 35  0.058333605 0.038162859 0.085848411 0.112428208 0.083850680 0.14938280
#> 45  0.038426386 0.022879981 0.062888388 0.069750070 0.046646521 0.10317765
#> 55  0.021432108 0.011446513 0.039560089 0.037504206 0.022144657 0.06300306
#> 65  0.009341936 0.004452949 0.019444989 0.016196526 0.008357390 0.03098127
#> 75  0.002900235 0.001232958 0.006865364 0.004813645 0.002173955 0.01055622
#> all 0.066061937 0.046211884 0.089938716 0.122187713 0.097842813 0.15247074
#>             all      all_lo      all_up
#> 5   0.258461625 0.190515615 0.351867131
#> 15  0.179695903 0.146646728 0.221528993
#> 25  0.125271008 0.103982952 0.151551642
#> 35  0.085321209 0.064913699 0.111754703
#> 45  0.054272143 0.036484548 0.080522327
#> 55  0.029739716 0.017595310 0.050301611
#> 65  0.013048353 0.006737108 0.025101124
#> 75  0.003995321 0.001788337 0.008811119
#> all 0.094655007 0.079077809 0.113988108
#> 
#> $ear
#>            male     male_lo    male_up      female   female_lo  female_up
#> 5   0.180991413 0.122464921 0.26521423 0.280928669 0.208319078 0.37917933
#> 15  0.145671293 0.103712020 0.20037124 0.225743123 0.182005184 0.27817512
#> 25  0.114078178 0.080903922 0.15644111 0.177243542 0.145031848 0.21589577
#> 35  0.084801146 0.057078481 0.12117249 0.132892978 0.102336908 0.17162829
#> 45  0.057506510 0.035995649 0.08941087 0.091672552 0.063786725 0.13028683
#> 55  0.034070636 0.019410638 0.05825693 0.055614271 0.034629985 0.08819714
#> 65  0.016666232 0.008590505 0.03177567 0.027746277 0.015417861 0.04936936
#> 75  0.005878352 0.002685266 0.01247551 0.009852212 0.004853038 0.01972337
#> all 0.083709520 0.059624009 0.11489726 0.124576903 0.102115845 0.15191113
#>             all      all_lo     all_up
#> 5   0.230527318 0.174563760 0.30513382
#> 15  0.184994650 0.153053809 0.22427950
#> 25  0.145300155 0.121361676 0.17432024
#> 35  0.108838660 0.084640990 0.13973037
#> 45  0.074955600 0.052345410 0.10666407
#> 55  0.045241776 0.028286413 0.07204407
#> 65  0.022642128 0.012592451 0.04051401
#> 75  0.008132142 0.003999695 0.01641731
#> all 0.104612606 0.087724086 0.12591100
#> 
```
