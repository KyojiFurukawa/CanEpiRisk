#'Calculate the cumulative excess risk due to radiation exposure
#'
#'@param exposure a list object, exposure scenario (a list object, which contains age(s) at exposure (a single value or a vector), doseGy in Gy or Sv (a single value or a vector),
#'@param reference baseline rate and all cause mortality rate in the reference population (a list object, which contains data.frame objects named 'baseline' for baseline rates of the target endpoint and 'mortality' for all cause mortality rates in the reference population)
#'@param riskmodel risk model risk model (a list object, which contains two list objects for excess relative risk model (err) and excess absolute risk model (ear), each of which contains a vector of parameter values (para), a matrix of variance covariance matrix (var), and a function to compute the risk given a parameter vector, a dose value, an age at exposure, an attained age and sex.
#'@param option option for risk calculation (a list object, which contains maximum age to follow up (an integer value) )
#'
#'@return risk information(vector)
#'
#'@examples
#'    # The following examples use default data provided in the CanEpiRisk package
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
#'@export
YLL <- function( exposure, reference, riskmodel, option )
{
  # exposure=list( agex=5, doseGy=0.1, sex=1 ); riskmodel=LSS_mortality$allsolid$L; reference=list( baseline=mortality_Japan2018$allsolid, mortality=mortality_Japan2018$allcause ); option=list( mc_para=NULL, maxage=100, err_wgt=1, n_mcsamp=10000)

  mle_YLL  <- mc_YLL( exposure, reference, riskmodel, option=list(mle_only=T, maxage=option$maxage, err_wgt=option$err_wgt) )
  samp <- mc_YLL(  exposure, reference, riskmodel, option )
  c( mle=mle_YLL, mean=mean(samp), median=median(samp), ci=quantile(samp,c(0.025,0.975)) )
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

population_YLL <- function( dsGy, reference, riskmodel, agex=1:8*10-5, nmc=10000 ){
  # dsGy=0.1; reference=ref0; riskmodel=rm0
  lars0 <- mc_popYLL( dsGy, riskmodel, reference, agexs=agex, n_mcsamp=nmc )
  list( err=popLAR( lars0=lars0, wgt=c(1,0), agedist=reference$agedist, PER=1, agex=agex ) ,
        ear=popLAR( lars0=lars0, wgt=c(0,1), agedist=reference$agedist, PER=1, agex=agex ) )
}

