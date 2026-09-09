# Cancer mortality rates of WHO global regions

Cancer mortality rates of WHO global regions.

## Usage

``` r
Mortality
```

## Format

A list object of the five WHO global regions:

- `Aus-NZ Europe Northern America`:

  a list object of 18 items, each of which contains a data.frame of
  cancer rates (see Details)

- `Northern Africa - Western Asia`:

  a list object of 18 items, each of which contains a data.frame of
  cancer rates (see Details)

- `Latin America and Caribbean`:

  a list object of 18 items, each of which contains a data.frame of
  cancer rates (see Details)

- `Asia excl. Western Asia`:

  a list object of 18 items, each of which contains a data.frame of
  cancer rates (see Details)

- `Sub-Saharan Africa`:

  a list object of 18 items, each of which contains a data.frame of
  cancer rates (see Details)

## Details

The list object for each region contains the 18 site-specific cancer
mortality rates (`esophagus`, `stomach`, `colon`, `liver`, `pancreas`,
`lung`, `breast`, `prostate`, `bladder`, `brainCNS`, `thyroid`,
`all_leukaemia`, `all_cancer`, `allsolid-NMSC`, `allsolid`, `leukaemia`,
`allcause`, `survival`).

Each site-specific data.frame contains variables `age`, `male` and
`female`.

For `allcause` mortality data.frame contains age-specific person-year
data, `male_py` and `female_py`, in addition to `age`, `male` and
`female`.

## See also

[Incidence](https://kyojifurukawa.github.io/CanEpiRisk/reference/Incidence.md),
[plot_refdata](https://kyojifurukawa.github.io/CanEpiRisk/reference/plot_refdata.md)

## Examples

``` r
 names(Mortality)       # WHO global regions
#> [1] "Aus-NZ Europe Northern America" "Northern Africa - Western Asia"
#> [3] "Latin America and Caribbean"    "Asia excl. Western Asia"       
#> [5] "Sub-Saharan Africa"            
 names(Mortality[[1]])  # Sites for which baseline mortality rates are available
#>  [1] "esophagus"     "stomach"       "colon"         "liver"        
#>  [5] "pancreas"      "lung"          "breast"        "prostate"     
#>  [9] "bladder"       "brainCNS"      "thyroid"       "all_leukaemia"
#> [13] "all_cancer"    "allsolid-NMSC" "allsolid"      "leukaemia"    
#> [17] "allcause"      "survival"     

 # Example 1: All solid cancer mortality rates for Region-1
 head( Mortality[[1]]$allsolid )
#>   age         male       female
#> 1   1 3.993729e-06 3.302123e-06
#> 2   2 1.198119e-05 9.906370e-06
#> 3   3 1.996865e-05 1.651062e-05
#> 4   4 1.946148e-05 1.625391e-05
#> 5   5 1.895432e-05 1.599719e-05
#> 6   6 1.844715e-05 1.574048e-05
 tail( Mortality[[1]]$allsolid )
#>     age       male     female
#> 95   95 0.02189662 0.01180741
#> 96   96 0.02189662 0.01180741
#> 97   97 0.02189662 0.01180741
#> 98   98 0.02189662 0.01180741
#> 99   99 0.02189662 0.01180741
#> 100 100 0.02189662 0.01180741

 # Example 2: Leukaemia mortality rates for Region-3
 head( Mortality[[3]]$leukaemia )
#>   age         male       female
#> 1   1 4.066002e-06 3.390723e-06
#> 2   2 1.219801e-05 1.017217e-05
#> 3   3 2.033001e-05 1.695362e-05
#> 4   4 2.103463e-05 1.704913e-05
#> 5   5 2.173925e-05 1.714464e-05
#> 6   6 2.244387e-05 1.724015e-05
 tail( Mortality[[3]]$leukaemia )
#>     age         male       female
#> 95   95 0.0004455293 0.0002765286
#> 96   96 0.0004455293 0.0002765286
#> 97   97 0.0004455293 0.0002765286
#> 98   98 0.0004455293 0.0002765286
#> 99   99 0.0004455293 0.0002765286
#> 100 100 0.0004455293 0.0002765286

 # Example 3: A;;ll;-cause mortality rates for Region-5
 head( Mortality[[5]]$allcause )
#>   age        male      female  male_py female_py
#> 1   1 0.056160846 0.048467260 18761.24  18253.13
#> 2   2 0.009220674 0.007515785 18127.85  17684.38
#> 3   3 0.007183781 0.006129537 17634.31  17232.95
#> 4   4 0.005848985 0.005303171 17194.27  16823.33
#> 5   5 0.004864717 0.004698722 16794.60  16448.30
#> 6   6 0.004105606 0.004193329 16422.67  16095.09
 tail( Mortality[[5]]$allcause )
#>     age      male    female male_py female_py
#> 95   95 0.3576011 0.2754637  9.8210   23.4550
#> 96   96 0.3680964 0.2902003  6.9710   16.8780
#> 97   97 0.3780983 0.3040167  4.9220   11.9500
#> 98   98 0.3884477 0.3176900  3.4625    8.3635
#> 99   99 0.3989252 0.3305587  2.4190    5.8265
#> 100 100 0.4149599 0.3458980  1.6845    4.0590

 # Example 4: plotting lung cancer mortality rates
 plot_refdata( dat=Mortality, outcome="lung", leg_pos=c(0.27,0.95) )

```
