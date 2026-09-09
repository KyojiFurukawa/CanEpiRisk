# population_LAR: Population-averaged Lifetime Attributable Risk

Compute **population-averaged lifetime attributable risk (LAR)** due to
radiation exposure, aggregating over a population age distribution and
sex-specific baseline rates. The function combines user-specified excess
risk models (ERR/EAR) with site-specific baseline incidence/mortality
and all-cause mortality to produce age-category summaries and overall
totals, with uncertainty from Monte Carlo sampling.

## Usage

``` r
population_LAR(
  dsGy,
  reference,
  riskmodel,
  agex = 1:8 * 10 - 5,
  PER = 10000,
  nmc = 10000
)
```

## Arguments

- dsGy:

  Numeric scalar. Radiation **dose** in Gy (or Sv if your model is
  parameterized accordingly). Must be nonnegative.

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
      models), *or* `ci` — confidence bounds for one-parameter models;

    - `f` — function of the form `f(beta, data, lag)` returning
      age-specific excess risk given parameters and data (dose, ages,
      sex).

  See package examples (e.g., `LSS_mortality$allsolid$L`,
  `LSS_incidence$leukaemia$LQ`).

- agex:

  Numeric vector of **ages at exposure** (midpoints for grouped
  categories). Default is `1:8 * 10 - 5` (i.e., 5, 15, …, 75)
  representing 0–10, 10–20, …, 70–80.

- PER:

  Integer scaling factor for reporting results *per* `PER` persons.
  Default `10000` (i.e., cases per 10,000 persons).

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

Values are LAR **per `PER` persons**.

## Details

Let `agedist` denote the population age distribution (by single-year or
grouped ages) and `baseline`/`mortality` denote site-specific baseline
rates and all-cause mortality for the same region/population and sex
coding. For each age-at-exposure category in `agex`, `population_LAR()`
evaluates LAR under the supplied risk model(s) and then averages across
the age structure, returning sex-specific estimates and an `all`-sex
average. Uncertainty is obtained by drawing model parameters either from
a variance–covariance matrix (`var`) or, in one-parameter models, from
95% confidence bounds (`ci`) when provided.

Results are reported **per `PER` persons** (default: per 10,000). The
default `agex = 1:8 * 10 - 5` corresponds to midpoint ages 5, 15, …, 75,
i.e., exposure categories 0–10, 10–20, …, 70–80 years, which should
match the categories of `agedist`.

## Units & Alignment

Doses must be in Gy (or consistent with your model). Ensure `baseline`
and `mortality` share the same age grid (typically 1–100) and sex
coding, and that `agedist` corresponds to the same population.

## See also

[`CER`](https://kyojifurukawa.github.io/CanEpiRisk/reference/CER.md),
[`plot_agedist`](https://kyojifurukawa.github.io/CanEpiRisk/reference/plot_agedist.md),
[`population_YLL`](https://kyojifurukawa.github.io/CanEpiRisk/reference/population_YLL.md),
[`YLL`](https://kyojifurukawa.github.io/CanEpiRisk/reference/YLL.md)

## Examples

``` r
set.seed(100)
# Default package data: LSS-derived models and WHO-like regional rates
# Example 1: All solid cancer mortality, Region 1, dose = 0.1 Gy
ref1 <- list(
  baseline  = Mortality[[1]]$allsolid,   # site-specific baseline mortality
  mortality = Mortality[[1]]$allcause,   # all-cause mortality
  agedist   = agedist_rgn[[1]]           # population age distribution
)
mod1 <- LSS_mortality$allsolid$L         # linear ERR model (LSS)
out1 <- population_LAR(dsGy = 0.1, reference = ref1, riskmodel = mod1)
out1$err   # per 10,000 persons by default
#>          male   male_lo   male_up     female  female_lo female_up        all
#> 5   130.37052 83.684260 203.04478 221.901462 145.741428 338.04836 175.529433
#> 15   92.10210 63.056550 131.64627 156.473982 114.794414 214.09618 123.449312
#> 25   65.42755 45.114384  91.52862 110.051941  85.437062 143.21353  87.395346
#> 35   46.41875 30.530566  68.26803  76.956627  57.227143 103.48037  61.729352
#> 45   32.54654 19.564183  52.98648  51.984849  35.094862  77.12590  42.443120
#> 55   21.22069 11.404601  38.98402  33.086807  19.583245  55.58879  27.321900
#> 65   12.22755  5.847979  25.35129  18.738471   9.682757  35.86489  15.754229
#> 75    5.77911  2.450553  13.67047   8.481243   3.835635  18.60684   7.335593
#> all  52.81522 36.466146  73.98011  83.455376  64.177472 109.88703  68.485470
#>         all_lo    all_up
#> 5   120.542859 259.25148
#> 15   95.165421 164.04028
#> 25   70.012219 110.64680
#> 35   46.722220  81.66375
#> 45   28.712903  62.47583
#> 55   16.167233  46.05064
#> 65    8.134999  30.36719
#> 75    3.280844  16.19603
#> all  54.678094  87.25822

# Example 2: Leukaemia incidence, Region 4, dose = 0.1 Gy, LQ model
ref2 <- list(
  baseline  = Incidence[[4]]$leukaemia,
  mortality = Mortality[[4]]$allcause,
  agedist   = agedist_rgn[[4]]
)
mod2 <- LSS_incidence$leukaemia$LQ
out2 <- population_LAR(dsGy = 0.1, reference = ref2, riskmodel = mod2, PER = 10000, nmc = 10000)
out2$err
#>          male    male_lo  male_up    female  female_lo female_up       all
#> 5   16.887526 0.09510633 71.64065 12.803411 0.07139758  52.82134 14.948397
#> 15   8.101229 0.03876649 20.11160  6.519586 0.03078176  15.82578  7.351688
#> 25   6.343635 0.03221517 13.38227  5.493388 0.02795455  11.57930  5.934571
#> 35   6.222154 0.02525119 12.92559  5.466720 0.02237215  11.48464  5.852195
#> 45   6.521093 0.02523791 14.57651  5.490056 0.02127415  12.34995  6.011522
#> 55   6.806484 0.02594877 17.28338  5.358925 0.02039229  13.58449  6.084538
#> 65   6.629825 0.02477937 19.66861  4.865219 0.01810706  14.22901  5.720060
#> 75   5.538629 0.02033106 19.23206  3.817298 0.01393391  12.94007  4.608737
#> all  8.807671 0.04953417 22.81574  6.954832 0.03980346  17.48186  7.897506
#>         all_lo   all_up
#> 5   0.08377032 62.67576
#> 15  0.03496200 18.09759
#> 25  0.03018724 12.51534
#> 35  0.02384939 12.21162
#> 45  0.02328657 13.48848
#> 55  0.02317223 15.43496
#> 65  0.02134263 16.86451
#> 75  0.01687491 15.83090
#> all 0.04477201 20.19862
```
