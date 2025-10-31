#' @title Calculate Years of Life Lost (YLL) due to Radiation Exposure
#'
#' @description
#' Computes the expected **Years of Life Lost (YLL)** attributable to radiation exposure.
#' This function integrates site-specific baseline mortality and all-cause mortality rates
#' with user-specified excess risk models (ERR/EAR) to estimate the expected reduction
#' in life expectancy due to radiation-attributable cancer mortality.
#'
#' The function supports both single-age and protracted exposure scenarios and allows
#' for uncertainty quantification via Monte Carlo sampling.
#'
#' @param exposure A **list** specifying the exposure scenario:
#'   \itemize{
#'     \item \code{agex} — age(s) at exposure (scalar or numeric vector).
#'     \item \code{doseGy} — dose(s) in Gray (Gy), scalar or vector of same length as \code{agex}.
#'     \item \code{sex} — sex indicator (1 = male, 2 = female).
#'   }
#'
#' @param reference A **list** specifying baseline and all-cause mortality information:
#'   \itemize{
#'     \item \code{baseline} — data.frame of site-specific mortality (or incidence) rates
#'       per person-year, including columns \code{age}, \code{male}, and \code{female}.
#'     \item \code{mortality} — data.frame of all-cause mortality rates (per person-year)
#'       with same structure and age range as \code{baseline}.
#'   }
#'   These should correspond to the same population or region.
#'
#' @param riskmodel A **list** defining the excess risk model, typically of the form:
#'   \itemize{
#'     \item \code{err} and/or \code{ear} — sublists each containing:
#'       \itemize{
#'         \item \code{para} — vector of model parameters.
#'         \item \code{var} — variance–covariance matrix of parameters (for multi-parameter models).
#'         \item \code{ci} — 95% confidence bounds for a single-parameter model (alternative to \code{var}).
#'         \item \code{f(beta, data, lag)} — function returning age-specific excess risks.
#'       }
#'   }
#'
#' @param option A **list** of optional settings:
#'   \itemize{
#'     \item \code{maxage} — maximum attained age for accumulation (default: 100).
#'     \item \code{err_wgt} — weight between ERR (1) and EAR (0) transfer models.
#'     \item \code{n_mcsamp} — number of Monte Carlo draws for uncertainty estimation (default: 10,000).
#'     \item \code{alpha} — significance level for confidence intervals (default: 0.05).
#'   }
#'
#' @details
#' The YLL estimate integrates the excess mortality over age, weighted by remaining life expectancy,
#' producing an interpretable measure of life expectancy loss attributable to radiation.
#' The computation uses region- and sex-specific mortality data provided in the package or supplied by the user.
#'
#' Typical usage involves pairing \code{YLL()} with \code{CER()} or \code{population_YLL()} for
#' complementary assessments of probability-based and life loss–based risk metrics.
#'
#' @return
#' A named numeric vector with point and interval summaries of cumulative excess risk,
#' typically including:
#' \itemize{
#'   \item \code{mle}: point estimate,
#'   \item \code{mean}, \code{median}: Monte Carlo summaries,
#'   \item \code{ci_lo}, \code{ci_up}: confidence interval.
#' }
#' Values are per person; multiply by \code{1e4} or \code{1e5} to report per 10,000 or 100,000.
#'
#'
#' @examples
#' set.seed(100)
#' ## Example 1: All-solid cancer mortality (Region 1, female, 0.1 Gy at age 15)
#' exp1 <- list(agex = 15, doseGy = 0.1, sex = 2)
#' ref1 <- list(
#'   baseline  = Mortality[[1]]$allsolid,
#'   mortality = Mortality[[1]]$allcause
#' )
#' mod1 <- LSS_mortality$allsolid$L
#' opt1 <- list(maxage = 100, err_wgt = 1, n_mcsamp = 10000)
#'
#' YLL(exposure = exp1, reference = ref1, riskmodel = mod1, option = opt1)
#'
#' ## Example 2: Leukaemia incidence (Region 4), male, 100 mGy evenly across ages 30–45
#' exp2 <- list(agex = 30:44 + 0.5, doseGy = rep(0.1 / 15, 15), sex = 1)
#' ref2 <- list(
#'   baseline  = Incidence[[4]]$leukaemia,
#'   mortality = Mortality[[4]]$allcause
#' )
#' mod2 <- LSS_incidence$leukaemia$LQ
#' opt2 <- list(maxage = 60, err_wgt = 0, n_mcsamp = 10000)
#'
#' YLL(exposure = exp2, reference = ref2, riskmodel = mod2, option = opt2)
#'
#' @seealso
#' \code{\link{CER}}, \code{\link{population_YLL}}, \code{\link{population_LAR}}
#'@importFrom MASS mvrnorm
#'@export
YLL <- function( exposure, reference, riskmodel, option )
{
  # exposure=list( agex=5, doseGy=0.1, sex=1 ); riskmodel=LSS_mortality$allsolid$L; reference=list( baseline=mortality_Japan2018$allsolid, mortality=mortality_Japan2018$allcause ); option=list( mc_para=NULL, maxage=100, err_wgt=1, n_mcsamp=10000)
  if( is.null(option$alpha) ) option$alpha <- 0.05  # alpha error (required to compute the confidence interval )
  mle_YLL  <- mc_YLL( exposure, reference, riskmodel, option=list(mle_only=T, maxage=option$maxage, err_wgt=option$err_wgt) )
  samp <- mc_YLL(  exposure, reference, riskmodel, option )
  ci <- quantile(samp,c(option$alpha/2,1-option$alpha/2))
  c( mle=mle_YLL, mean=mean(samp), median=median(samp), ci_lo=ci[1], ci_up=ci[2] )
}

#'Calculating population-averaged years of life lost
#'
#' @description
#' Compute **population-averaged years of life lost (YLL)** due to radiation
#' exposure, aggregating over a population age distribution and sex-specific baseline
#' rates. The function combines user-specified excess risk models (ERR/EAR) with
#' site-specific baseline incidence/mortality and all-cause mortality to produce
#' age-category summaries and overall totals, with uncertainty from Monte Carlo
#' sampling.
#'
#' @details
#' Let \code{agedist} denote the population age distribution (grouped
#' ages) and \code{baseline}/\code{mortality} denote site-specific baseline rates and
#' all-cause mortality for the same region/population and sex coding. For each
#' age-at-exposure category in \code{agex}, \code{population_YLL()} evaluates YLL under
#' the supplied risk model(s) and then averages across the age structure, returning
#' sex-specific estimates and an \code{all}-sex average. Uncertainty is obtained by
#' drawing model parameters either from a variance–covariance matrix (\code{var}) or,
#' in one-parameter models, from 95% confidence bounds (\code{ci}) when provided.
#'
#' The default
#' \code{agex = 1:8 * 10 - 5} corresponds to midpoint ages 5, 15, …, 75, i.e.,
#' exposure categories 0–10, 10–20, …, 70–80 years, which should match the categories of \code{agedist}.
#'
#' @param dsGy Numeric scalar. Radiation **dose** in Gy (or Sv). Must be nonnegative.
#'
#' @param reference List describing the target **reference population**, with components:
#'   \itemize{
#'     \item \code{baseline} — data frame of site-specific baseline **incidence or mortality**
#'       rates (per person-year) over ages 1–100, with columns \code{age}, \code{male}, \code{female}.
#'     \item \code{mortality} — data frame of **all-cause mortality** (same structure/age grid).
#'     \item \code{agedist} — data frame or vector with the **population age distribution** used
#'       to average risks across ages (e.g., by single year or grouped ages).
#'   }
#'   All three components must refer to the **same region/population** and share consistent
#'   sex coding.
#'
#' @param riskmodel List defining the radiation **risk model**. Typically contains
#'   sublists for ERR and/or EAR:
#'   \itemize{
#'     \item \code{err} / \code{ear} — each is a list with:
#'       \itemize{
#'         \item \code{para} — numeric vector of model parameters;
#'         \item \code{var} — parameter variance–covariance matrix (for multi-parameter models), \emph{or}
#'               \code{ci} — confidence interval for one-parameter models;
#'         \item \code{f} — function of the form \code{f(beta, data, lag)} returning age-specific
#'               excess risk given parameters and data (dose, ages, sex).
#'       }
#'   }
#'   See package examples (e.g., \code{LSS_mortality$allsolid$L}, \code{LSS_incidence$leukaemia$LQ}).
#'
#' @param agex Numeric vector of **ages at exposure** (midpoints for grouped categories).
#'   Default is \code{1:8 * 10 - 5} (i.e., 5, 15, …, 75) representing 0–10, 10–20, …, 70–80.
#'
#' @param nmc Integer Monte Carlo sample size for uncertainty propagation. Default \code{10000}.
#'
#' @return
#' A **list** with components for each transfer model present (typically \code{$err} and \code{$ear}).
#' Each component is a data frame with one row per exposure age category (rows named by
#' \code{agex}) plus an \code{all} row, and columns:
#' \itemize{
#'   \item \code{male}, \code{male_lo}, \code{male_up}
#'   \item \code{female}, \code{female_lo}, \code{female_up}
#'   \item \code{all}, \code{all_lo}, \code{all_up}
#' }
#' Values are YLL.
#'
#' @section Units & Alignment:
#' Doses must be in Gy (or consistent with your model). Ensure \code{baseline} and
#' \code{mortality} share the same age grid (typically 1–100) and sex coding, and that
#' \code{agedist} corresponds to the same population.
#'
#' @examples
#' set.seed(100)
#' # The following examples use default data provided in CanEpiRisk package
#' # for riskmodels (LSS_mortality) derived from Life Span Study
#' #     baseline rates and age distribution for WHO riskmodels (Mortality, agedist_rgn)
#' # Example: allsolid mortality, Region-1, exposed to 0.1 Gy, followed up to age 100, LSS linear ERR
#' ref1 <- list(  baseline=Mortality[[1]]$allsolid,     # baseline rates
#' mortality=Mortality[[1]]$allcause,     # allcause mortality
#' agedist=agedist_rgn[[1]] )           # age distribution
#' mod1 <- LSS_mortality$allsolid$L                     # risk model
#' population_YLL( dsGy=0.1, reference=ref1, riskmodel=mod1 )
#'
#' @seealso
#' \code{\link{YLL}}, \code{\link{plot_agedist}}, \code{\link{population_CER}}, \code{\link{CER}}
#'@importFrom MASS mvrnorm
#'@export
population_YLL <- function( dsGy, reference, riskmodel, agex=1:8*10-5, nmc=10000 ){
  # dsGy=0.1; reference=ref0; riskmodel=rm0
  lars0 <- mc_popYLL( dsGy, riskmodel, reference, agexs=agex, n_mcsamp=nmc )
  list( err=popLAR( lars0=lars0, wgt=c(1,0), agedist=reference$agedist, PER=1, agex=agex ) ,
        ear=popLAR( lars0=lars0, wgt=c(0,1), agedist=reference$agedist, PER=1, agex=agex ) )
}


mc_YLL <- function( exposure, reference, riskmodel, option )
{
  # option=list( mc_para=NULL, maxage=100, err_wgt=1, n_mcsamp=10000)
  # option=list(mle_only=T, maxage=option$maxage, err_wgt=option$err_wgt)
  ages <- reference$baseline$age
  nexp <- length(exposure$agex)
  sexlab <- c("male", "female")[exposure$sex]
  mrate <- reference$mortality[[sexlab]]         # all cause
  brate <- reference$baseline[[sexlab]]          # baseline cancer rate
  survp <- c(1, exp(-cumsum(mrate))[-length(mrate)])

  ok <- min(exposure$agex) < ages & ages <= option$maxage

  if( !is.null(option$mle_only) ){
    if( option$mle_only ){
      a <- sapply(1:nexp, function(i) riskmodel$err$f( riskmodel$err$para,
                                                       data.frame(dose=exposure$doseGy[i], agex=exposure$agex[i], age=ages-0.5, sex=exposure$sex) ) )
      err <- apply(a, 1, sum)
      brate_d <- brate * (1+err)
      mrate_d <- brate_d + (mrate - brate)
      survp_d <- c(1, exp(-cumsum(mrate_d))[-length(mrate_d)])
      mle_yll_err <- sum( (survp-survp_d)[ok]/survp[min(exposure$agex) + 1] )

      a <- sapply(1:nexp, function(i) riskmodel$ear$f( riskmodel$ear$para,
                                                       data.frame(dose=exposure$doseGy[i], agex=exposure$agex[i], age=ages-0.5, sex=exposure$sex) ) )
      ear <- apply(a, 1, sum)

      brate_d <- brate + ear
      mrate_d <- brate_d + (mrate - brate)
      survp_d <- c(1, exp(-cumsum(mrate_d))[-length(mrate_d)])
      mle_yll_ear <- sum( (survp-survp_d)[ok]/survp[min(exposure$agex) + 1] )

      return( option$err_wgt * mle_yll_err + (1-option$err_wgt) * mle_yll_ear )
    }
  }

  n_mcsamp <- option$n_mcsamp
  mc_para  <- option$mc_para
  if( !is.null(mc_para$err) ) n_mcsamp <- nrow( mc_para$err )

  YLLs_err <- YLLs_ear <- matrix( 0, nrow=sum(ok), ncol=n_mcsamp )

  if( option$err_wgt > 0 ){
    if( is.null(mc_para$err) ) mc_para$err <- MASS::mvrnorm( n=n_mcsamp, mu=riskmodel$err$para, Sigma=riskmodel$err$var )
    sim_err <- apply( mc_para$err, 1, function(bet) { # bet=mc_para$err[1,]
      a <- sapply(1:nexp, function(i) riskmodel$err$f( bet, data.frame(dose=exposure$doseGy[i], agex=exposure$agex[i], age=ages-0.5, sex=exposure$sex)))
      apply(a, 1, sum)    })
    brate_d <- brate * (1+sim_err)
    mrate_d <- brate_d + (mrate - brate)   #  mrate + brate*sim_err
    survp_d <- apply( mrate_d, 2, function(x) c(1, exp(-cumsum(x))[-length(x)]) )
    YLLs_err <- apply( survp_d, 2, function(spd) { (survp-spd)[ok]/survp[min(exposure$agex) + 1] } )
  }

  if( option$err_wgt < 1 ){
    if( is.null(mc_para$ear) ) mc_para$ear <- MASS::mvrnorm( n=n_mcsamp, mu=riskmodel$ear$para, Sigma=riskmodel$ear$var )
    sim_ear <- apply(mc_para$ear, 1, function(bet) {
      a <- sapply(1:nexp, function(i) riskmodel$ear$f(bet, data.frame(dose=exposure$doseGy[i], agex=exposure$agex[i], age=ages-0.5, sex=exposure$sex)))
      apply(a, 1, sum)     })
    brate_d <- brate + sim_ear
    mrate_d <- brate_d + (mrate - brate)
    survp_d <- apply( mrate_d, 2, function(x) c(1, exp(-cumsum(x))[-length(x)]) )
    YLLs_ear <- apply(survp_d, 2, function(spd) { (survp-spd)[ok]/survp[min(exposure$agex) + 1] });
  }
  res_byage <- option$err_wgt * YLLs_err + (1-option$err_wgt) * YLLs_ear
  res_total <- apply(res_byage, 2, sum)
  res_total
}


mc_popYLL <- function( ds, riskmodel, reference, agexs, n_mcsamp=0 ){
  # ds=0.1; riskmodel=rm0; reference=ref0; n_mcsamp=B; agexs=AGEX

  if( n_mcsamp > 0 ){
    mc_paras_err <- MASS::mvrnorm( n=n_mcsamp, mu=riskmodel$err$para, Sigma=riskmodel$err$var )
    mc_paras_ear <- MASS::mvrnorm( n=n_mcsamp, mu=riskmodel$ear$para, Sigma=riskmodel$ear$var )
    mc_para0 <- list( err=mc_paras_err, ear=mc_paras_ear )
  } else {
    mc_para0 <- list( err=matrix(riskmodel$err$para,byrow=T,nrow=1), ear=matrix(riskmodel$ear$para,byrow=T,nrow=1) )
  }
  option0   <- list( mc_para=mc_para0, maxage=100, err_wgt=1, n_mcsamp=NULL)   # ERR

  ylls_m_err <- sapply( agexs, function(ax) mc_YLL( list(agex=ax,doseGy=ds,sex=1), reference, riskmodel, option=option0 ) )
  ylls_f_err <- sapply( agexs, function(ax) mc_YLL( list(agex=ax,doseGy=ds,sex=2), reference, riskmodel, option=option0 ) )

  option0$err_wgt <- 0  # EAR
  ylls_m_ear <- sapply( agexs, function(ax) mc_YLL( list(agex=ax,doseGy=ds,sex=1), reference, riskmodel, option=option0 ) )
  ylls_f_ear <- sapply( agexs, function(ax) mc_YLL( list(agex=ax,doseGy=ds,sex=2), reference, riskmodel, option=option0 ) )

  list( err=list(male=ylls_m_err, female=ylls_f_err),
        ear=list(male=ylls_m_ear, female=ylls_f_ear)  )
}


