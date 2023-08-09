#'Calculate the cumulated excess risk due to radiation exposure
#'
#'@param exposure exposure scenario (a list object, which contains age(s) at exposure (a single value or a vector), doseGy in Gy or Sv (a single value or a vector),
#'@param doseGy dose in Gy or Sv (a single value or a vector)
#'@param sex sex 1:male 2:female
#'@param maxage maximum age to follow up (an integer value)
#'@param riskmodel risk model (a list object, which contains two list objects for excess relative risk model (err) and excess absolute risk model (ear), each of which contains a vector of parameter values (para), a matrix of variance covariance matrix (var), and a function to compute the risk given a parameter vector, a dose value, an age at exposure, an attained age and sex.
#'@param wgt weights for ERR vs EAR transfer (e.g. c(1,0) (default), which indicates ERR transfer)
#'@param baseline age- and sex-specific baseline rate (data.frame with columns for age, male, female); default: all solid cancer mortality rate of Japan in 2018
#'@param mortality age- and sex-specific all cause mortality rate (data.frame with columns for age, male, female); default: all cause mortality rate of Japan in 2018
#'@param alpha significance level for the confidence interval (default value = 0.05, which corresponds to 95\% confidence interval)
#'@param n_mcsamp number of Monte Carlo sample size (default:10000)
#'@param seed random number seed (a single value; default is no seed)
#'
#'@return risk information(vector)
#'
#'@examples
#'    # The following examples use default data provided in the SUMRAY package
#'    # for riskmodel (LSS R14 all solid cancer model),
#'    #     baseline (all solid cancer mortality rates in Japan 2018) and
#'    #     mortality  (all cause mortality rates in Japan 2018)
#'
#'    # Cumulated excess risk for male exposed to 0.1 Gy at age 10
#'    #  followed up to age 90 with err transfer
#'    CER( agex=10, doseGy=0.1, sex=1, maxage=90 )
#'
#'    # Cumulated excess risk for female exposed to 0.01 Gy/year at ages 10-19
#'    #  followed up to age 100 with 7:3 weights for ERR and EAR transfers
#'    CER( agex=10:19, doseGy=rep(0.01,10), sex=2, maxage=100, wgt=c(.7,.3) )
#'
#'@importFrom MASS mvrnorm
#'#'@export
mc_CER <- function( exposure, reference, riskmodel, option )
{
  #   exposure=exposure0; reference=reference0; riskmodel=riskmodel0; option=option0
  ages <- reference$baseline$age
  nexp <- length(exposure$agex)
  sexlab <- c("male", "female")[exposure$sex]
  mrate <- reference$mortality[[sexlab]]
  brate <- reference$baseline[[sexlab]]
  survp <- c(1, exp(-cumsum(mrate))[-length(mrate)])
  n_mcsamp <- option$n_mcsamp
  mc_para <- option$mc_para
  if( !is.null(mc_para$err) ) n_mcsamp <- nrow( mc_para$err )

  ok <- min(exposure$agex) < ages & ages <= option$maxage
  CERs_err <- CERs_ear <- matrix( 0, nrow=sum(ok), ncol=n_mcsamp )

  if( option$wgt[1] > 0 ){
    if( is.null(mc_para$err) ) mc_para$err <- MASS::mvrnorm( n=n_mcsamp, mu=riskmodel$err$para, Sigma=riskmodel$err$var )
    sim_err <- apply( mc_para$err, 1, function(bet) { # bet=mc_para$err[1,]
      a <- sapply(1:nexp, function(i) riskmodel$err$f( bet, data.frame(dose=exposure$doseGy[i], agex=exposure$agex[i], age=ages, sex=exposure$sex)))
      apply(a, 1, sum)
    })
    CERs_err <- apply(sim_err, 2, function(err) { (brate * err * survp)[ok]/survp[min(exposure$agex) + 1] })
  }
  if( option$wgt[2] > 0 ){
    if( is.null(mc_para$ear) ) mc_para$ear <- MASS::mvrnorm( n=n_mcsamp, mu=riskmodel$ear$para, Sigma=riskmodel$ear$var )
    sim_ear <- apply(mc_para$ear, 1, function(bet) {
      a <- sapply(1:nexp, function(i) riskmodel$ear$f(bet, data.frame(dose=exposure$doseGy[i], agex=exposure$agex[i], age=ages, sex=exposure$sex)))
      apply(a, 1, sum)
    })
    CERs_ear <- apply(sim_ear, 2, function(ear) { (ear * survp)[ok]/survp[min(exposure$agex) + 1] })
  }

  res_byage <- option$wgt[1] * CERs_err + option$wgt[2] * CERs_ear
  res_total <- apply(res_byage, 2, sum)
  res_total
}


CER <- function( exposure, reference, riskmodel, option )
{
  samp_CER <- mc_CER(  exposure, reference, riskmodel, option )
  c( mean=mean(samp_CER), median=median(samp_CER), ci=quantile(samp_CER,c(0.025,0.975)) )
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

  option0   <- list( mc_para=mc_para0, maxage=100, wgt=c(1,0), n_mcsamp=NULL)   # ERR

  n_agexs <- length( agexs )
  lars_m_err <- sapply( agexs, function(ax) mc_CER( list(agex=ax,doseGy=ds,sex=1), reference, riskmodel, option=option0 ) )
  lars_f_err <- sapply( agexs, function(ax) mc_CER( list(agex=ax,doseGy=ds,sex=2), reference, riskmodel, option=option0 ) )

  option0$wgt <- c(0,1)  # EAR
  lars_m_ear <- sapply( agexs, function(ax) mc_CER( list(agex=ax,doseGy=ds,sex=1), reference, riskmodel, option=option0 ) )
  lars_f_ear <- sapply( agexs, function(ax) mc_CER( list(agex=ax,doseGy=ds,sex=2), reference, riskmodel, option=option0 ) )

  list( err=list(male=lars_m_err, female=lars_f_err),
        ear=list(male=lars_m_ear, female=lars_f_ear)  )
}
#--------------------------------------#
# poppulation calculation              #
#--------------------------------------#
popLAR <- function( lars0, wgt, agedist, PER=100 ){     #    lars0 <- lars_inci_leukaemia_LQ_100; wgt=c(1,0)   # 表にまとめる
  wlars0 <- list( male=NULL, female=NULL )
  wlars0 <- list(     male = wgt[1]*lars0$err$male   + wgt[2]*lars0$ear$male,
                      female = wgt[1]*lars0$err$female + wgt[2]*lars0$ear$female )

  tab1 <- cbind( agex=AGEX, t( apply( wlars0$male,   2, quantile, prob=c(0.5, 0.025, 0.975) ) )*PER,     # convert to %
                 t( apply( wlars0$female, 2, quantile, prob=c(0.5, 0.025, 0.975) ) )*PER )

  ( lar_mf <- rbind(   male=quantile( apply( t(wlars0$male)   * agedist$p_mf[,"male"],   2, sum ), prob=c(0.5,0.025,0.975) ),
                       female=quantile( apply( t(wlars0$female) * agedist$p_mf[,"female"], 2, sum ), prob=c(0.5,0.025,0.975) ) ) * PER )
  ( lar_age <- t( apply( t(wlars0$male) * agedist$p_age[,"male"] + t(wlars0$female) * agedist$p_age[,"female"], 1, quantile, prob=c(0.5, 0.025, 0.975) )* PER) )
  ( lar_all <- quantile( apply( t(wlars0$male)*agedist$p[,"male"] + t(wlars0$female)*agedist$p[,"female"], 2, sum), prob=c(0.5, 0.025, 0.975) )* PER )

  ret <- rbind( cbind( tab1[,2:7], lar_age ),
                c( lar_mf[1,], lar_mf[2,], lar_all ) )
  rownames(ret) <- c( as.character(AGEX), "all" )
  colnames(ret) <- c( "male", "male_lo", "male_up", "female", "female_lo", "female_up", "all", "all_lo", "all_up" )
  ret
}


population_LAR <- function( dsGy, reference, riskmodel, AGEX=1:8*10-5  ){    # dsGy=0.1; reference=ref0; riskmodel=rm0
  lars0 <- mc_popLAR( dsGy, riskmodel, reference, agexs=AGEX, n_mcsamp=B )
  list( err=popLAR( lars0=lars0, wgt=c(1,0), agedist=reference$agedist, PER=100 ) ,
        ear=popLAR( lars0=lars0, wgt=c(0,1), agedist=reference$agedist, PER=100 ) )
}

