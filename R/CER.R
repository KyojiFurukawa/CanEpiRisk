#'Calculate the cumulative excess risk due to radiation exposure
#'
#'@param exposure a list object that specifies the exposure scenario, which contains 'agex' (a single value or a vector for age(s) at exposure), 'doseGy' (a single value or a vector of dose(s) in Gy), and 'sex' (1 or 2 for male or female).
#'@param reference a list object that specifies the baseline information of the reference population, which contains data.frame objects named 'baseline' for baseline rates of the target endpoint and 'mortality' for all cause mortality rates.
#'@param riskmodel a list object that specifies the risk model, which contains two list objects named 'err' for excess relative rate model and 'ear' for excess absolute rate model, each of which contains a vector 'para' for model parameter estimates, a matrix 'var' for the variance covariance matrix, and a function 'f' to compute the excess risk given a parameter vector and exposure information (e.g., dose, age at exposure, sex, attained age).
#'@param option a list object that specifies optional settings for risk calculation, which contains an integer value 'maxage' for the maximum age to follow up, a value 'err_wgt' for the weight for risk transfer (1=err, 0=ear), an integer value 'n_mcsamp' for the number of Monte Carlo samples, and alpha for the significance level (default=0.05).
#'
#'@return information of calculated risk (vector)
#'
#'@examples
#'    # The following examples use default data provided in the CanEpiRisk package
#'    # allsolid mortality, Region-1, male, 0.1Gy at age 5, LSS linear ERR
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
#'#'@export
population_LAR <- function( dsGy, reference, riskmodel, agex=1:8*10-5, PER=100, nmc=10000 ){    # dsGy=0.1; reference=ref0; riskmodel=rm0
  lars0 <- mc_popLAR( dsGy, riskmodel, reference, agexs=agex, n_mcsamp=nmc )
  list( err=popLAR( lars0=lars0, wgt=c(1,0), agedist=reference$agedist, PER=PER, agex=agex ) ,
        ear=popLAR( lars0=lars0, wgt=c(0,1), agedist=reference$agedist, PER=PER, agex=agex ) )
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

  # S_t <- survp/survp[min(exposure$agex) + 1]
  # S_t[S_t>1] <- 0
  # send.table.to.excel( data.frame( ages-1, ages, ages-0.5, brate, mrate, survp, err, brate*err, S_t,  brate*err*S_t) )

  #  par( mfrow=c(1, 5 ) )
  #  plot( ages, brate, type="l", xlab="age" )
  #  plot( ages[ages>5], err[ages>5], type="l", xlim=c(0,100), xlab="age" )
  #  plot( ages[ages>5], (brate*err)[ages>5], type="l", xlim=c(0,100), xlab="age" )
  #  plot( ages[ages>5], S_t[ages>5], type="l", xlim=c(0,100), xlab="age" )
  #  plot( ages[ages>5], (brate*err*S_t)[ages>5], type="l", xlim=c(0,100), xlab="age" )

  #  sum( brate*err*S_t )

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



popLAR <- function( lars0, wgt, agedist, PER=100, agex ){     #    lars0 <- lars_inci_leukaemia_LQ_100; wgt=c(1,0)
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






