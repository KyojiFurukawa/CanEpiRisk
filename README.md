# CanEpiRisk

# R Library for Calculation of Lifetime Risks due to Radiation Exposure

Version: 1.0-0 Date: 2025-10-28

Author & Maintainer: Kyoji Furukawa 

URL: <https://github.com/KyojiFurukawa/CanEpiRisk>

## Overview

Quantifying cancer risks associated with ionizing radiation exposure is a central concern in radiation epidemiology and public health, particularly when estimating lifetime risks across populations and exposure scenarios. `CanEpiRisk` is a comprehensive R package designed to facilitate radiation-associated lifetime risk calculations. The package provides an integrated framework for importing exposure and demographic data, transferring specified risk models—such as **excess relative risk (ERR)** and **excess absolute risk (EAR)** and computing key risk measures, including **cumulative excess risk (CER)** and **years of life lost (YLL)**. By offering a unified interface for evaluating risks under various exposure scenarios and across diverse population settings, **CanEpiRisk** streamlines the workflow for radiation risk assessment, enhances analytical reproducibility, and enables transparent comparison of model-based estimates.

---


## Installation

```         
# install.packages("devtools")
devtools::install_github("KyojiFurukawa/CanEpiRisk")
library(CanEpiRisk)
```


## Quick Start

[Using or specifying reference data](fbi.pdf).

[Using or specifying risk models](https://github.com/KyojiFurukawa/CanEpiRisk/vignettes/CanEpiRisk_Models.html)

[Specifying exposure scenarios](fbi.pdf).

## Quick examples

```         
# Example 1: allsolid mortality, Region-1, female, 0.1Gy at age 15, followed up to age 100, LSS linear ERR
exp1 <- list( agex=5, doseGy=0.1, sex=2 )   # exposure scenario
ref1 <- list( baseline=Mortality[[1]]$allsolid,        # baseline rates
             mortality=Mortality[[1]]$allcause )       # all-cause mortality
mod1 <- LSS_mortality$allsolid$L                       # risk model
opt1 <- list( maxage=100, err_wgt=1, n_mcsamp=10000 )  # option
CER(  exposure=exp1, reference=ref1, riskmodel=mod1, option=opt1 ) * 10000 # cases per 10,000

# Example 2: leukaemia incidence, Region-4, male, 6.7(100/15)mGy at ages 30-45, followed up to age 60, LSS LQ EAR
exp2 <- list( agex=30:44+0.5, doseGy=rep(0.1/15,15), sex=1 )
ref2 <- list( baseline=Incidence[[4]]$leukaemia,       # baseline rates
             mortality=Mortality[[4]]$allcause )       # all-cause mortality rates
mod2 <- LSS_incidence$leukaemia$LQ                     # risk model
opt2 <- list( maxage=60, err_wgt=0, n_mcsamp=10000)    # option
CER(  exposure=exp2, reference=ref2, riskmodel=mod2, option=opt2 ) * 10000 # cases per 10,000
```





