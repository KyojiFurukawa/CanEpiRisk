#' Calculating Cumulative Excess Risks (CER)
#'
#' @title CER: Calculating Cumulative Excess Risks
#'
#' @description
#' Compute the **cumulative excess risk (CER)** attributable to radiation exposure for an
#' individual scenario or cohort, integrating excess risks over age up to a specified
#' maximum while accounting for competing all-cause mortality and model uncertainty.
#' This is the core engine for lifetime risk estimation in **CanEpiRisk**.
#' :contentReference[oaicite:0]{index=0}
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
#'     \code{var} (variance–covariance matrix) \emph{or} \code{ci} (length-2 95\% CI for
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
#'   \item \code{ci.2.5\%}, \code{ci.97.5\%}: 95% interval bounds.
#' }
#' Values are per person; multiply by \code{1e4} or \code{1e5} to report per 10,000 or 100,000.
#'
#' @section Units & Alignment:
#' Doses must be in Gy. Ensure \code{baseline} and \code{mortality} tables are from the same
#' population/region and aligned on age (1–100) and sex coding. :contentReference[oaicite:1]{index=1}
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
#'@seealso \link{population_LAR}, \link{YLL}
#'@importFrom MASS mvrnorm
#'@export
CER <- function( exposure, reference, riskmodel, option )
{
  if( is.null(option$alpha) ) option$alpha <- 0.05  # alpha error (required to compute the confidence interval )
  mle_CER  <- mc_CER(  exposure, reference, riskmodel, option=list(mle_only=T, maxage=option$maxage, err_wgt=option$err_wgt ) )
  if( is.null( riskmodel$err$var )  ){   # in case of a riskmodel with a single err/ear parameter
    rm <- riskmodel
    rm$err$para <- rm$err$ci[1]; rm$ear$para <- rm$ear$ci[1]
    loci_CER  <- mc_CER(  exposure, reference, rm, option=list(mle_only=T, maxage=option$maxage, err_wgt=option$err_wgt ) )
    rm$err$para <- rm$err$ci[2]; rm$ear$para <- rm$ear$ci[2]
    upci_CER  <- mc_CER(  exposure, reference, rm, option=list(mle_only=T, maxage=option$maxage, err_wgt=option$err_wgt ) )
    res <- c( mle=mle_CER, mean=mle_CER, median=mle_CER, ci=c(loci_CER, upci_CER) )
  } else {
    samp_CER <- mc_CER(  exposure, reference, riskmodel, option )
    res <- c( mle=mle_CER, mean=mean(samp_CER), median=median(samp_CER), ci=quantile(samp_CER,c(option$alpha/2,1-option$alpha/2)) )
  }
  res
}

#'Calculating population-averaged lifetime attributable risks
#'@description Calculate the population-averaged lifetime attributable risk due to radiation exposure.
#'
#'@param dsGy radiation dose in Gy or Sv (a single value).
#'@param reference baseline rate, all cause mortality rate and age distribution in the reference population (a list object, which contains data.frame objects named \code{baseline} for baseline rates of the target endpoint, 'mortality' for all cause mortality rates and 'agedist' for age distribution in the reference population).
#'@param riskmodel risk model risk model (a list object, which contains two list objects for excess relative risk model (\code{err}) and excess absolute risk model (\code{ear}), each of which contains a vector of parameter values (para), a matrix of variance covariance matrix (var), and a function to compute the risk given a parameter vector, a dose value, an age at exposure, an attained age and sex.
#'@param agex a vector of ages at exposure, which represent age categories (default values: 5, 15, 25, ..., 75 to represent age categories 0-10, 10-20, ..., 70-80) option for risk calculation (a list object, which contains maximum age to follow up (an integer value)
#'@param PER  an integer (default value: 10000 to show the estimated risk as cases per 10000)
#'@param mmc  an integer for the Monte Carlo sample size (default: 10000)
#'
#'@return estimated risk information (list)
#'
#'@examples
#'  # The following examples use default data provided in CanEpiRisk package
#'  # for riskmodels (LSS_mortality and LSS_incidence) derived from Life Span Study
#'  # and baseline mortality and incidence rates for WHO global regions (Mortality and Incidence).
#'
#'  # Example 1: allsolid mortality, Region-1, exposed to 0.1 Gy, followed up to age 100, LSS linear ERR
#'  ref1 <- list(  baseline=Mortality[[1]]$allsolid,     # baseline rates
#'                mortality=Mortality[[1]]$allcause,     # allcause mortality
#'                  agedist=agedist_rgn[[1]] )           # age distribution
#'  mod1 <- LSS_mortality$allsolid$L                     # risk model
#'  population_LAR( dsGy=0.1, reference=ref1, riskmodel=mod1 )    # CER cases per 10,000
#'
#'  # Example 2: leukaemia incidence, Region-4, exposed to 0.1 Gy, followed up to age 100, LSS LQ ERR
#'  ref2 <- list(  baseline=Incidence[[4]]$leukaemia,    # baseline rates
#'                mortality=Mortality[[4]]$allcause,     # all-cause mortality
#'                  agedist=agedist_rgn[[4]] )           # age distribution
#'  mod2 <- LSS_incidence$leukaemia$LQ                   # risk model
#'  population_LAR( dsGy=0.1, reference=ref2, riskmodel=mod2 )    # CER cases per 10,000
#'
#'@seealso \link{CER}, \link{population_YLL}
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





mc_CER <- function( exposure, reference, riskmodel, option )
{
  #   exposure=exposure0; reference=reference0; riskmodel=riskmodel0; option=option0
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
  # ds=0.1; riskmodel=rm0; reference=ref0; n_mcsamp=B; agexs=AGEX

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













