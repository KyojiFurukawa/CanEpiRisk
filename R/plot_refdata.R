#'Plotting reference data
#'@description plot_refdata Plots the reference data of cancer incidence or mortality rates for the five WHO global regions.
#'
#'@param dat Mortality or Incidence, which contains age- and sex-specific cancer incidence and mortality rates by cancer site and region.
#'@param outcome a character string that specifies the cancer site for which the rates are plotted.
#'@param title a character string that  specifies the title of the plot.
#'@param ymax a value that specifies the maximum value for y-axis.
#'@param leg_pos a vector that specifies the legend x-y position ((1,1) for the top-right; default="none" for no-legend).
#'@param PER an integer value for the rate denominator (default=10^5)
#'@param x_lab x-axis label (default="age (years)")
#'@param y_lab y-axis label (default="cases per 100,000 PY")
#'@param lsz   legend size (default=9)
#'
#'@return a ggplot object
#'
#'@details
#' The parameter outcome can be chosen from one of the following character strings for both Mortality and Incidence:
#'   "esophagus", "stomach", "colon", "liver", "pancreas", "lung", "breast", "prostate", "bladder"
#'   "brainCNS", "thyroid", "all_leukaemia", "all_cancer", "allsolid-NMSC", "allsolid",
#'   "leukaemia" (leukaemia excluding CLL).
#'
#' For MOrtality, "allcause" and "survival" can be additionally chosen.
#'
#'
#'@examples
#'  # The following examples use default data provided in CanEpiRisk package
#'  # for riskmodels (LSS_mortality and LSS_incidence) derived from Life Span Study
#'  # and baseline mortality and incidence rates for WHO global regions (Mortality and Incidence)
#'
#'  # Example 1: All solid cancer mortality rates
#'  plot_refdata( dat=Mortality, outcome="allsolid", title="All solid cancers mortality", leg_pos=c(0.27,0.95), y_lab="deaths per 100,000 PY" )
#'
#'  # Example 2: Leukaemia incidence rates
#'  plot_refdata( dat=Incidence, outcome="leukaemia", title="Leukaemia incidence (excluding CLL)", leg_pos=c(0.27,0.95), y_lab="cases per 100,000 PY" )
#'
#'@importFrom ggplot2 ggplot aes geom_line xlab ylab facet_grid margin ggtitle scale_color_manual element_rect element_text theme_bw theme
#'@export
plot_refdata <- function( dat, outcome, title="", ymax=NULL, leg_pos="none", PER=100000,
                          x_lab="age (years)", y_lab="cases per 100,000 PY", lsz=9 ){
  a <- NULL
  for( i in 1:length(dat) ){ b <- dat[[i]][[outcome]]
  a <- rbind( a, data.frame( outcome=outcome, region=factor( names(Mortality)[i], levels=names(Mortality)[i] ),
                             sex=factor(rep( c("male","female"), each=nrow(b) ), levels=c("male","female")),
                             age=rep( b$age, 2 ),
                             rate=c(b$male,b$female)*PER  ) ) }

  if( is.null(ymax) ) ymax <- max(a$rate)
  g <- ggplot( a, aes( x=age, y=rate ) ) + geom_line( ggplot2::aes(col=region) ) +
    xlab(x_lab) + ylab(y_lab)  + facet_grid( ~sex, scale="free") + ggtitle(title) +
    scale_color_manual( values=c("black","red","green","blue","orange") )  + theme_bw() +
    theme(legend.position=leg_pos, legend.justification=c("right","top"), legend.box.just="right", legend.margin=margin(1,1,1,1),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.direction = "vertical", legend.text = element_text(size=lsz), legend.title = element_text(size=0)  )
  g
}

#'Plotting risk models
#'@description plot_riskmodel() plots the excess relative and absolute rates of cancer incidence or mortality rates based on the risk models derived from the Life Span Study for a given dose. The scenario is currently fixed for ages at exposure of 10, 30 and 50 years. Age specific risks are plotted for both ERR and EAR separately for both male and female.
#'@param rm a list object which contains the risk model information specified in the format described in Details.
#'@param doseGy a value that specifies the radiation dose in Gy or Sv (default=0.1).
#'@param maxage maximum age in years (default=100)
#'@param ymax a value that specifies the maximum value for y-axis.
#'@param leg_pos a vector that specifies the legend x-y position ((1,1) for the top-right; default="none" for no-legend).
#'@param EAR_PER an integer value for the excess absolute rate denominator (default=10^4).
#'@param add a vector of values to be added to female's ERR and EAR values (to avoid completely overlapping sex-specific values when the risk is not modified by sex).
#'@param sex a character string which specifies whether the plot contains "both" sexes, "male" or "female" only.
#'@param x_lab x-axis label (default="age (years)" ).
#'@param y_lab y-axis label (default="excess rate").
#'@param title a character string that specifies the title of the plot (default="", i.e., no title).
#'
#'@return a ggplot object
#'
#'@details
#' The risk model can be specified by the list object provided in CanEpiRisk package; for example
#'     LSS_mortality$allsolid$L  for LSS all solid cancer mortality linear model,
#'     LSS_incidence$allsolid$LQ for LSS all solid cancer incidence linear-quadratic model,
#'     LSS_mortality$leukaemia$LQ for LSS leukaemia mortality linear-quadratic model,
#'     LSS_incidence$thyroid$L for LSS thyroid cancer incidence linear model.
#'
#'
#'@examples
#'  # The following examples use default data provided in CanEpiRisk package
#'  # for riskmodels (LSS_mortality and LSS_incidence) derived from Life Span Study
#'  # and baseline mortality and incidence rates for WHO global regions (Mortality and Incidence)
#'
#'  # Example 1: LSS all solid cancer mortality risk model
#'  plot_riskmodel( rm=LSS_mortality$allsolid$L, title="LSS all solid cancer mortality, Linear",  leg_pos=c(0.4, 0.95) )
#'
#'  # Example 2: LSS Leukaemia incidence risk model
#'  plot_riskmodel( rm=LSS_incidence$leukaemia$LQ, title="LSS leukaemia incidence", ymax=c(1.5, .3), add=c(0.01,0) )
#'
#'@importFrom ggplot2 ggplot aes geom_line xlab ylab facet_wrap scale_y_continuous ggtitle geom_point theme_bw margin  scale_color_manual element_rect element_text theme_bw theme
#'@export
plot_riskmodel <- function( rm, doseGy=0.1, maxage=100, ymax=c(1.5, 30), leg_pos="none", EAR_PER=10^4, add=c(0,0), sex="both",
                            x_lab="age (years)", y_lab="excess rate", title="" ){
  # rm=LSS_incidence$esophagus$L; doseGy=0.1; maxage=100; modelid=""; ymax=NULL;EAR_PER=1;add=c(0,0)
  a1m <- data.frame( age=15:maxage, sex="male", age_at_exp=10 ); a1f <- data.frame( age=15:maxage, sex="female", age_at_exp=10 )
  a2m <- data.frame( age=35:maxage, sex="male", age_at_exp=30 ); a2f <- data.frame( age=35:maxage, sex="female", age_at_exp=30)
  a3m <- data.frame( age=55:maxage, sex="male", age_at_exp=50 ); a3f <- data.frame( age=55:maxage, sex="female", age_at_exp=50 )
  aaa <- data.frame( rbind( a1m, a1f, a2m, a2f, a3m, a3f, a1m, a1f, a2m, a2f, a3m, a3f ) )
  aaa$age_at_exp=factor( aaa$age_at_exp )
  aaa$sex <- factor( aaa$sex, levels=c("male","female") )
  aaa$transfer=factor( rep(c("ERR","EAR"), each=nrow(aaa)/2 ) , levels=c("ERR","EAR" ) )
  aaa$risk=c( rm$err$f( beta=rm$err$para, data.frame( dose=doseGy, agex=10, age=15:maxage, sex=1) ),
              rm$err$f( beta=rm$err$para, data.frame( dose=doseGy, agex=10, age=15:maxage, sex=2) )+ add[1],
              rm$err$f( beta=rm$err$para, data.frame( dose=doseGy, agex=30, age=35:maxage, sex=1) ),
              rm$err$f( beta=rm$err$para, data.frame( dose=doseGy, agex=30, age=35:maxage, sex=2) )+ add[1],
              rm$err$f( beta=rm$err$para, data.frame( dose=doseGy, agex=50, age=55:maxage, sex=1) ),
              rm$err$f( beta=rm$err$para, data.frame( dose=doseGy, agex=50, age=55:maxage, sex=2) )+ add[1],
              rm$ear$f( beta=rm$ear$para, data.frame( dose=doseGy, agex=10, age=15:maxage, sex=1) )*EAR_PER,
              rm$ear$f( beta=rm$ear$para, data.frame( dose=doseGy, agex=10, age=15:maxage, sex=2) )*EAR_PER + add[2],
              rm$ear$f( beta=rm$ear$para, data.frame( dose=doseGy, agex=30, age=35:maxage, sex=1) )*EAR_PER,
              rm$ear$f( beta=rm$ear$para, data.frame( dose=doseGy, agex=30, age=35:maxage, sex=2) )*EAR_PER + add[2],
              rm$ear$f( beta=rm$ear$para, data.frame( dose=doseGy, agex=50, age=55:maxage, sex=1) )*EAR_PER,
              rm$ear$f( beta=rm$ear$para, data.frame( dose=doseGy, agex=50, age=55:maxage, sex=2) )*EAR_PER + add[2])
  bbb <-  expand.grid( sex=unique(aaa$sex), transfer=unique(aaa$transfer) )

  if( is.null(ymax) ) ymax <- c( max( aaa$risk[ aaa$transfer=="ERR"] ) ,  max( aaa$risk[ aaa$transfer=="EAR" ] ) )

  if( sex=="male" ) aaa$risk[ aaa$sex=="female" ] <- NA
  if( sex=="female" ) aaa$risk[ aaa$sex=="male" ] <- NA

  bbb$risk <- rep(ymax, each=2)
  bbb$age <- 20
  bbb$age_at_exp <- aaa$age_at_exp[1]

  aaa$risk[ aaa$transfer=="ERR" & aaa$risk > ymax[1] ] <- NA
  aaa$risk[ aaa$transfer=="EAR" & aaa$risk > ymax[2] ] <- NA
  g <- ggplot2::ggplot( aaa, aes( x=age, y=risk, color=sex, linetype=age_at_exp ) ) + facet_wrap(.~transfer, scale="free") + scale_y_continuous(limits = c(0, NA)) +
    ylab(y_lab) + xlab(x_lab) + ggtitle(title) +
    geom_line()  + scale_color_manual(values=c("blue","red") ) + theme_bw() + geom_point(data=bbb, alpha=0) +
    theme(legend.position = leg_pos, legend.justification=c("right","top"), legend.box.just="right", legend.margin=margin(1,1,1,1),
          legend.direction = "vertical", legend.text = element_text(size=8), legend.title = element_text(size=8)  )
  g
  #  + scale_color_manual("sex",c("male"="blue","female"="red") )
}

#'Plotting age distributions
#'@description plot_agedist() plots the age distribution of the WHO global regions.
#'@param region a single value or a vector which specifies the WHO gloval region(s) (default=1:5; 1="Aus-NZ Europe Northern America", 2="Northern Africa - Western Asia", 3="Latin America and Caribbean", 4="Asia excl. Western Asia", 5="Sub-Saharan Africa")
#'@param agedist a list object which contains information for the region-, sex- and age-specific distribution (default is "agedist_rgn" object of CanEpiRisk)
#'@param rgn_labs a vector of character strings which specifies the labels of the region names shown in the plot.
#'
#'@return a ggplot object
#'
#'
#'@examples
#'  # The following examples use default data provided in CanEpiRisk package
#'  # for age distribution for WHO riskmodels (agedist_rgn)
#'
#'  # Example: age distributions for Regions 1 and 5
#'  plot_agedist( regions=c(1,5) )
#'
#'@importFrom ggplot2 ggplot aes geom_col xlab ylab xlim facet_grid scale_fill_manual scale_x_continuous
#'@export
plot_agedist <- function( regions=1:5, agedist=agedist_rgn, rgn_labs=names(agedist_rgn) ){
  pdat <- NULL
  for( rgn in regions ){
    a <- agedist[[rgn]][[1]]
    b <- data.frame( Sex=factor(rep(c("male","female"), each=8), levels=c("male","female")), Age=rownames(a),
                     population_percentage=unlist(a)*100, Region=rgn_labs[rgn] )
    b$population_percentage[9:16] <- -b$population_percentage[9:16]
    pdat <- rbind( pdat, b )
  }
  pdat$Region <- factor(pdat$Region, levels=unique(pdat$Region) )
  ( pop_range <- range( pdat$population_percentage ) )
  pop_range_seq <- seq(from=pop_range[1], to=pop_range[2], length = 6)
  pop_range_seq
  ( pop_range_breaks <- pretty(pop_range, n = 7) )

  ggplot( pdat, aes( x=population_percentage, y=Age, fill=Sex ) ) + geom_col(width=0.5) + xlab("Population (%)") +
    facet_grid( ~Region, scale="free") + scale_fill_manual(values=c("blue", "red")) + xlim(c(-15,15)) +
    scale_x_continuous(breaks=pop_range_breaks,
                       labels=abs(pop_range_breaks))
}







