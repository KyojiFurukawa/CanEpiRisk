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
#'   \item \code{ci.2.5\%}, \code{ci.97.5\%}: 95% interval bounds.
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
#' ## Example 2: Protracted exposure (EAR, leukaemia incidence)
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
  c( mle=mle_YLL, mean=mean(samp), median=median(samp), ci=quantile(samp,c(option$alpha/2,1-option$alpha/2)) )
}

#'Calculating population-averaged years of life lost
#'@description Calculate the population-averaged years of life lost due to radiation exposure.
#'
#'@param dsGy radiation dose in Gy or Sv (a single value).
#'@param reference baseline rate, all cause mortality rate and age distribution in the reference population (a list object, which contains data.frame objects named 'baseline' for baseline rates of the target endpoint, 'mortality' for all cause mortality rates and 'agedist' for age distribution in the reference population).
#'@param riskmodel risk model risk model (a list object, which contains two list objects for excess relative risk model (err) and excess absolute risk model (ear), each of which contains a vector of parameter values (para), a matrix of variance covariance matrix (var), and a function to compute the risk given a parameter vector, a dose value, an age at exposure, an attained age and sex.
#'@param agex a vector of ages at exposure, which represent age categories (default values: 5, 15, 25, ..., 75 to represent age categories 0-10, 10-20, ..., 70-80) option for risk calculation (a list object, which contains maximum age to follow up (an integer value)
#'@param mmc  an integer for the Monte Carlo sample size (default: 10000)
#'
#'@return estimated risk information (list)
#'
#'@examples
#'  set.seed(100)
#'  # The following examples use default data provided in CanEpiRisk package
#'  # for riskmodels (LSS_mortality) derived from Life Span Study
#'  #     baseline rates and age distribution for WHO riskmodels (Mortality, agedist_rgn)
#'
#'  # Example: allsolid mortality, Region-1, exposed to 0.1 Gy, followed up to age 100, LSS linear ERR
#'  ref1 <- list(  baseline=Mortality[[1]]$allsolid,     # baseline rates
#'                mortality=Mortality[[1]]$allcause,     # allcause mortality
#'                  agedist=agedist_rgn[[1]] )           # age distribution
#'  mod1 <- LSS_mortality$allsolid$L                     # risk model
#'  population_YLL( dsGy=0.1, reference=ref1, riskmodel=mod1 )    # YLL
#'
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


