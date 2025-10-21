
plot_refdata <- function( dat, outcome, title="", ymax=NULL, leg_pos="none", PER=100000,
                          x_lab="age (years)", y_lab="cases per 100,000 PY", lsz=9 ){
  a <- NULL
  for( i in 1:length(dat) ){ b <- dat[[i]][[outcome]]
  a <- rbind( a, data.frame( outcome=outcome, region=factor( names(Mortality)[i], levels=names(Mortality)[i] ),
                             sex=factor(rep( c("male","female"), each=nrow(b) ), levels=c("male","female")),
                             age=rep( b$age, 2 ),
                             rate=c(b$male,b$female)*PER  ) ) }

  if( is.null(ymax) ) ymax <- max(a$rate)
  g <- ggplot2::ggplot( a, aes( x=age, y=rate ) ) + geom_line( ggplot2::aes(col=region) ) +
    xlab(x_lab) + ylab(y_lab)  + facet_grid( ~sex, scale="free") + ggtitle(title) +
    scale_color_manual( values=c("black","red","green","blue","orange") )  + theme_bw() +
    theme(legend.position=leg_pos, legend.justification=c("right","top"), legend.box.just="right", legend.margin=margin(1,1,1,1),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.direction = "vertical", legend.text = element_text(size=lsz), legend.title = element_text(size=0)  )
  g
}

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
