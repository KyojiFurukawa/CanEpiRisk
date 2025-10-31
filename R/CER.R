#' Calculating Cumulative Excess Risks (CER)
#'
#' @title CER: Calculating Cumulative Excess Risks
#'
#' @description
#' Compute the **cumulative excess risk (CER)** attributable to radiation exposure for an
#' individual scenario or cohort, integrating excess risks over age up to a specified
#' maximum while accounting for competing all-cause mortality and model uncertainty.
#' This is the core engine for lifetime risk estimation in **CanEpiRisk**.
#'
#' @details
#' The function integrates age-specific excess risks implied by the supplied
#' ERR/EAR model(s) across the attained-age grid (typically ages 1–100). When
#' `n_mcsamp > 1`, parameter uncertainty is propagated via Monte Carlo sampling
#' using either the model variance–covariance matrix (`var`) or a 95% CI (`ci`)
#' for 1-parameter models. For protracted exposures (vectorized `agex`/`doseGy`),
#' contributions are summed across exposure segments. Baseline site-specific rates
#' and all-cause mortality must correspond to the same population/region and share
#' the same age grid and sex coding.
#'
#' @param exposure list. Exposure scenario with components:
#' \itemize{
#'   \item \code{agex}: age(s) at exposure (scalar or vector).
#'   \item \code{doseGy}: dose(s) in gray (Gy); same length as \code{agex} if vectorized.
#'   \item \code{sex}: sex indicator (\code{1} = male, \code{2} = female).
#' }
#'
#' @param reference list. Baseline reference data for the same population/region with:
#' \itemize{
#'   \item \code{baseline}: data frame of site-specific baseline rates (incidence or mortality),
#'         with columns \code{age}, \code{male}, \code{female} on ages 1:100.
#'   \item \code{mortality}: data frame of all-cause mortality (same columns/age grid).
#' }
#'
#' @param riskmodel list. Radiation risk model definition with sublists for
#' \emph{excess relative risk} (ERR) and \emph{excess absolute risk} (EAR), e.g.:
#' \itemize{
#'   \item \code{err}/\code{ear}: each a list containing
#'     \code{para} (numeric parameter vector),
#'     \code{var} (variance–covariance matrix) \emph{or} \code{ci} (confidence interval for
#'     1-parameter models), and
#'     \code{f} (function of the form \code{f(beta, data, lag)} returning age-specific excess risk).
#' }
#'
#' @param option list. Optional settings:
#' \itemize{
#'   \item \code{maxage}: maximum attained age for accumulation (e.g., \code{100}).
#'   \item \code{err_wgt}: weight to blend ERR vs EAR (\code{1} = pure ERR; \code{0} = pure EAR;
#'         intermediate values allowed).
#'   \item \code{n_mcsamp}: Monte Carlo sample size for uncertainty propagation (e.g., \code{10000}).
#'   \item \code{alpha}: significance level for interval estimation (default \code{0.05}).
#' }
#'
#' @return
#' A named numeric vector with point and interval summaries of cumulative excess risk,
#' typically including:
#' \itemize{
#'   \item \code{mle}: point estimate,
#'   \item \code{mean}, \code{median}: Monte Carlo summaries,
#'   \item \code{ci_lo}, \code{ci_up}: confidence intervals.
#' }
#' Values are per person; multiply by \code{1e4} or \code{1e5} to report per 10,000 or 100,000.
#'
#' @section Units & Alignment:
#' Doses must be in Gy. Ensure \code{baseline} and \code{mortality} tables are from the same
#' population/region and aligned on age (1–100) and sex coding.
#'
#' @examples
#' set.seed(100)
#' # Example 1: All-solid cancer mortality (Region 1), female, 0.1 Gy at age 15
#' exp1 <- list(agex = 15, doseGy = 0.1, sex = 2)
#' ref1 <- list(
#'   baseline  = Mortality[[1]]$allsolid,
#'   mortality = Mortality[[1]]$allcause
#' )
#' mod1 <- LSS_mortality$allsolid$L
#' opt1 <- list(maxage = 100, err_wgt = 1, n_mcsamp = 10000)
#'
#' CER(exposure = exp1, reference = ref1, riskmodel = mod1, option = opt1) * 10000
#'
#' # Example 2: Leukaemia incidence (Region 4), male, 100 mGy evenly across ages 30–45
#' exp2 <- list(agex = 30:44 + 0.5, doseGy = rep(0.1/15, 15), sex = 1)
#' ref2 <- list(
#'   baseline  = Incidence[[4]]$leukaemia,
#'   mortality = Mortality[[4]]$allcause
#' )
#' mod2 <- LSS_incidence$leukaemia$LQ
#' opt2 <- list(maxage = 60, err_wgt = 0, n_mcsamp = 10000)
#'
#' CER(exposure = exp2, reference = ref2, riskmodel = mod2, option = opt2) * 10000
#'
#'
#' @seealso
#' \code{\link{population_LAR}}, \code{\link{population_YLL}}, \code{\link{YLL}}
#'@importFrom MASS mvrnorm
#'@export
CER <- function( exposure, reference, riskmodel, option )
{
  if( is.null(option$alpha) ) option$alpha <- 0.05  # alpha error (required for specifying the level of confidence interval)
  mle_CER  <- mc_CER(  exposure, reference, riskmodel, option=list(mle_only=T, maxage=option$maxage, err_wgt=option$err_wgt ) )
  if( is.null( riskmodel$err$var )  ){  # If riskmodel has only a single err/ear parameter with a confidence interval...
    rm <- riskmodel
    rm$err$para <- rm$err$ci[1]; rm$ear$para <- rm$ear$ci[1]
    loci_CER  <- mc_CER(  exposure, reference, rm, option=list(mle_only=T, maxage=option$maxage, err_wgt=option$err_wgt ) )
    rm$err$para <- rm$err$ci[2]; rm$ear$para <- rm$ear$ci[2]
    upci_CER  <- mc_CER(  exposure, reference, rm, option=list(mle_only=T, maxage=option$maxage, err_wgt=option$err_wgt ) )
    res <- c( mle=mle_CER, mean=mle_CER, median=mle_CER, ci_lo=loci_CER, ci_up=upci_CER )
  } else {
    samp_CER <- mc_CER(  exposure, reference, riskmodel, option )
    ci <- c( quantile(samp_CER,c(option$alpha/2,1-option$alpha/2) ) )
    res <- c( mle=mle_CER, mean=mean(samp_CER), median=median(samp_CER), ci_lo=ci[1], ci_up=ci[2] )
  }
  res
}

#' @title population_LAR: Population-averaged Lifetime Attributable Risk
#'
#' @description
#' Compute **population-averaged lifetime attributable risk (LAR)** due to radiation
#' exposure, aggregating over a population age distribution and sex-specific baseline
#' rates. The function combines user-specified excess risk models (ERR/EAR) with
#' site-specific baseline incidence/mortality and all-cause mortality to produce
#' age-category summaries and overall totals, with uncertainty from Monte Carlo
#' sampling.
#'
#' @details
#' Let \code{agedist} denote the population age distribution (by single-year or grouped
#' ages) and \code{baseline}/\code{mortality} denote site-specific baseline rates and
#' all-cause mortality for the same region/population and sex coding. For each
#' age-at-exposure category in \code{agex}, \code{population_LAR()} evaluates LAR under
#' the supplied risk model(s) and then averages across the age structure, returning
#' sex-specific estimates and an \code{all}-sex average. Uncertainty is obtained by
#' drawing model parameters either from a variance–covariance matrix (\code{var}) or,
#' in one-parameter models, from 95% confidence bounds (\code{ci}) when provided.
#'
#' Results are reported **per \code{PER} persons** (default: per 10,000). The default
#' \code{agex = 1:8 * 10 - 5} corresponds to midpoint ages 5, 15, …, 75, i.e.,
#' exposure categories 0–10, 10–20, …, 70–80 years, which should match the categories of \code{agedist}.
#'
#' @param dsGy Numeric scalar. Radiation **dose** in Gy (or Sv if your model is parameterized
#'   accordingly). Must be nonnegative.
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
#'               \code{ci} — confidence bounds for one-parameter models;
#'         \item \code{f} — function of the form \code{f(beta, data, lag)} returning age-specific
#'               excess risk given parameters and data (dose, ages, sex).
#'       }
#'   }
#'   See package examples (e.g., \code{LSS_mortality$allsolid$L}, \code{LSS_incidence$leukaemia$LQ}).
#'
#' @param agex Numeric vector of **ages at exposure** (midpoints for grouped categories).
#'   Default is \code{1:8 * 10 - 5} (i.e., 5, 15, …, 75) representing 0–10, 10–20, …, 70–80.
#'
#' @param PER Integer scaling factor for reporting results \emph{per} \code{PER} persons.
#'   Default \code{10000} (i.e., cases per 10,000 persons).
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
#' Values are LAR **per \code{PER} persons**.
#'
#' @section Units & Alignment:
#' Doses must be in Gy (or consistent with your model). Ensure \code{baseline} and
#' \code{mortality} share the same age grid (typically 1–100) and sex coding, and that
#' \code{agedist} corresponds to the same population.
#'
#'
#' @examples
#' set.seed(100)
#' # Default package data: LSS-derived models and WHO-like regional rates
#' # Example 1: All solid cancer mortality, Region 1, dose = 0.1 Gy
#' ref1 <- list(
#'   baseline  = Mortality[[1]]$allsolid,   # site-specific baseline mortality
#'   mortality = Mortality[[1]]$allcause,   # all-cause mortality
#'   agedist   = agedist_rgn[[1]]           # population age distribution
#' )
#' mod1 <- LSS_mortality$allsolid$L         # linear ERR model (LSS)
#' out1 <- population_LAR(dsGy = 0.1, reference = ref1, riskmodel = mod1)
#' out1$err   # per 10,000 persons by default
#'
#' # Example 2: Leukaemia incidence, Region 4, dose = 0.1 Gy, LQ model
#' ref2 <- list(
#'   baseline  = Incidence[[4]]$leukaemia,
#'   mortality = Mortality[[4]]$allcause,
#'   agedist   = agedist_rgn[[4]]
#' )
#' mod2 <- LSS_incidence$leukaemia$LQ
#' out2 <- population_LAR(dsGy = 0.1, reference = ref2, riskmodel = mod2, PER = 10000, nmc = 10000)
#' out2$err
#'
#' @seealso
#' \code{\link{CER}}, \code{\link{plot_agedist}}, \code{\link{population_YLL}}, \code{\link{YLL}}
#'@importFrom MASS mvrnorm
#'@export
population_LAR <- function( dsGy, reference, riskmodel, agex=1:8*10-5, PER=10000, nmc=10000 ){
  lars0 <- mc_popLAR( dsGy, riskmodel, reference, agexs=agex, n_mcsamp=nmc )
  list( err=popLAR( lars0=lars0, wgt=c(1,0), agedist=reference$agedist, PER=PER, agex=agex ) ,
        ear=popLAR( lars0=lars0, wgt=c(0,1), agedist=reference$agedist, PER=PER, agex=agex ) )
}


#'Calculating excess risks
#'@description Calculate the excess risk from a risk model under a specified exposure scenario.
#'
#'@param exposure a list object that specifies the exposure scenario, which contains \code{agex} (a single value or a vector for age(s) at exposure), 'doseGy' (a single value or a vector of dose(s) in Gy), and 'sex' (1 or 2 for male or female).
#'@param riskmodel a list object that specifies the risk model, which contains two list objects named \code{err} for excess relative rate model and 'ear' for excess absolute rate model, each of which contains a vector 'para' for model parameter estimates and a function 'f' to compute the excess risk given a parameter vector and exposure information (e.g., dose, age at exposure, sex, attained age).
#'@param option a list object that specifies optional settings for risk calculation, which contains an integer value 'maxage' for the maximum age to follow up and a value 'err_wgt' for the weight for risk transfer (1=err, 0=ear).
#'@param per an integer value for the risk denominator (default=1).
#'@return information of calculated excess risk (data.frame)
#'
#'@examples
#'  # The following examples use default data provided in CanEpiRisk package
#'  # for riskmodels (LSS_mortality and LSS_incidence) derived from Life Span Study
#'  # and baseline mortality and incidence rates for WHO global regions (Mortality and Incidence).
#'
#'  # Example 1: allsolid mortality, Region-1, female, 0.1Gy at age 15, followed up to age 100, LSS linear ERR
#'  exp1 <- list( agex=5, doseGy=0.1, sex=2 )   # exposure scenario
#'  ref1 <- list( baseline=Mortality[[1]]$allsolid,        # baseline rates
#'               mortality=Mortality[[1]]$allcause )       # all-cause mortality
#'  mod1 <- LSS_mortality$allsolid$L                       # risk model
#'  opt1 <- list( maxage=100, err_wgt=1, n_mcsamp=10000 )  # option
#'  CER(  exposure=exp1, reference=ref1, riskmodel=mod1, option=opt1 ) * 10000 # cases per 10,000
#'
#'
#'@seealso \link{LSS_mortality}, \link{LSS_incidence}
#'@export
Comp_Exrisk <- function( exposure, riskmodel, option, per=1 ){
  ages <- ceiling(exposure$agex[1]):option$maxage
  nexp <- length(exposure$agex)
  rm <- riskmodel$err
  if( option$err_wgt==0 ) rm <- riskmodel$ear
  a <- sapply(1:nexp, function(i)
    rm$f( rm$para, data.frame(dose=exposure$doseGy[i], agex=exposure$agex[i],
                              age=ages-0.5, sex=exposure$sex) ) )
  b <- apply( a, 1, sum )
  data.frame( age=ages, risk=b*per )
}

mc_CER <- function( exposure, reference, riskmodel, option ){
  ages <- reference$baseline$age
  nexp <- length(exposure$agex)
  sexlab <- c("male", "female")[exposure$sex]
  mrate <- reference$mortality[[sexlab]]
  brate <- reference$baseline[[sexlab]]
  survp <- c(1, exp(-cumsum(mrate))[-length(mrate)])

  ok <- min(exposure$agex) < ages & ages <= option$maxage

  if( !is.null(option$mle_only) ){
    if( option$mle_only ){
      a <- sapply(1:nexp, function(i) riskmodel$err$f( riskmodel$err$para,
                                                       data.frame(dose=exposure$doseGy[i], agex=exposure$agex[i], age=ages-0.5, sex=exposure$sex) ) )
      err <- apply(a, 1, sum)
      mle_cer_err <- sum( (brate * err * survp)[ok]/survp[min(exposure$agex) + 1] )

      a <- sapply(1:nexp, function(i) riskmodel$ear$f( riskmodel$ear$para,
                                                       data.frame(dose=exposure$doseGy[i], agex=exposure$agex[i], age=ages-0.5, sex=exposure$sex) ) )
      ear <- apply(a, 1, sum)
      mle_cer_ear <- sum( ( ear * survp)[ok]/survp[min(exposure$agex) + 1] )

      return( option$err_wgt * mle_cer_err + (1-option$err_wgt) * mle_cer_ear )

    }
  }

  n_mcsamp <- option$n_mcsamp
  mc_para <- option$mc_para
  if( !is.null(mc_para$err) ) n_mcsamp <- nrow( mc_para$err )


  CERs_err <- CERs_ear <- matrix( 0, nrow=sum(ok), ncol=n_mcsamp )

  if( option$err_wgt > 0 ){
    if( is.null(mc_para$err) ) mc_para$err <- MASS::mvrnorm( n=n_mcsamp, mu=riskmodel$err$para, Sigma=riskmodel$err$var )
    sim_err <- apply( mc_para$err, 1, function(bet) { # bet=mc_para$err[1,]
      a <- sapply(1:nexp, function(i) riskmodel$err$f( bet, data.frame(dose=exposure$doseGy[i], agex=exposure$agex[i], age=ages-0.5, sex=exposure$sex)))
      apply(a, 1, sum)
    })
    CERs_err <- apply(sim_err, 2, function(err) { (brate * err * survp)[ok]/survp[min(exposure$agex) + 1] })
  }
  if( option$err_wgt < 1 ){
    if( is.null(mc_para$ear) ) mc_para$ear <- MASS::mvrnorm( n=n_mcsamp, mu=riskmodel$ear$para, Sigma=riskmodel$ear$var )
    sim_ear <- apply(mc_para$ear, 1, function(bet) {
      a <- sapply(1:nexp, function(i) riskmodel$ear$f(bet, data.frame(dose=exposure$doseGy[i], agex=exposure$agex[i], age=ages-0.5, sex=exposure$sex)))
      apply(a, 1, sum)
    })
    CERs_ear <- apply(sim_ear, 2, function(ear) { (ear * survp)[ok]/survp[min(exposure$agex) + 1] })
  }

  res_byage <- option$err_wgt * CERs_err + (1-option$err_wgt) * CERs_ear
  res_total <- apply(res_byage, 2, sum)
  res_total
}



mc_popLAR <- function( ds, riskmodel, reference, agexs, n_mcsamp=0 ){
  if( n_mcsamp > 0 ){
    mc_paras_err <- MASS::mvrnorm( n=n_mcsamp, mu=riskmodel$err$para, Sigma=riskmodel$err$var )
    mc_paras_ear <- MASS::mvrnorm( n=n_mcsamp, mu=riskmodel$ear$para, Sigma=riskmodel$ear$var )
    mc_para0 <- list( err=mc_paras_err, ear=mc_paras_ear )
  } else {
    mc_para0 <- list( err=matrix(riskmodel$err$para,byrow=T,nrow=1), ear=matrix(riskmodel$ear$para,byrow=T,nrow=1) )
  }

  option0   <- list( mc_para=mc_para0, maxage=100, err_wgt=1, n_mcsamp=NULL)   # ERR

  n_agexs <- length( agexs )
  lars_m_err <- sapply( agexs, function(ax) mc_CER( list(agex=ax,doseGy=ds,sex=1), reference, riskmodel, option=option0 ) )
  lars_f_err <- sapply( agexs, function(ax) mc_CER( list(agex=ax,doseGy=ds,sex=2), reference, riskmodel, option=option0 ) )

  option0$err_wgt <- 0  # EAR
  lars_m_ear <- sapply( agexs, function(ax) mc_CER( list(agex=ax,doseGy=ds,sex=1), reference, riskmodel, option=option0 ) )
  lars_f_ear <- sapply( agexs, function(ax) mc_CER( list(agex=ax,doseGy=ds,sex=2), reference, riskmodel, option=option0 ) )

  list( err=list(male=lars_m_err, female=lars_f_err),
        ear=list(male=lars_m_ear, female=lars_f_ear)  )
}



popLAR <- function( lars0, wgt, agedist, PER, agex ){     #    lars0 <- lars_inci_leukaemia_LQ_100; wgt=c(1,0)
  wlars0 <- list( male=NULL, female=NULL )
  wlars0 <- list(     male = wgt[1]*lars0$err$male   + wgt[2]*lars0$ear$male,
                      female = wgt[1]*lars0$err$female + wgt[2]*lars0$ear$female )

  tab1 <- cbind( agex=agex, t( apply( wlars0$male,   2, quantile, prob=c(0.5, 0.025, 0.975) ) )*PER,     # convert to %
                 t( apply( wlars0$female, 2, quantile, prob=c(0.5, 0.025, 0.975) ) )*PER )

  ( lar_mf <- rbind(   male=quantile( apply( t(wlars0$male)   * agedist$p_mf[,"male"],   2, sum ), prob=c(0.5,0.025,0.975) ),
                       female=quantile( apply( t(wlars0$female) * agedist$p_mf[,"female"], 2, sum ), prob=c(0.5,0.025,0.975) ) ) * PER )
  ( lar_age <- t( apply( t(wlars0$male) * agedist$p_age[,"male"] + t(wlars0$female) * agedist$p_age[,"female"], 1, quantile, prob=c(0.5, 0.025, 0.975) )* PER) )
  ( lar_all <- quantile( apply( t(wlars0$male)*agedist$p[,"male"] + t(wlars0$female)*agedist$p[,"female"], 2, sum), prob=c(0.5, 0.025, 0.975) )* PER )

  ret <- rbind( cbind( tab1[,2:7], lar_age ),
                c( lar_mf[1,], lar_mf[2,], lar_all ) )
  rownames(ret) <- c( as.character(agex), "all" )
  colnames(ret) <- c( "male", "male_lo", "male_up", "female", "female_lo", "female_up", "all", "all_lo", "all_up" )
  ret
}













