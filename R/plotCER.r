get_longtab <- function(obj0){   # obj0=Mortab$allsolid$L_100
  ( err <- obj0$err[1:8, ] )
  ( ear <- obj0$ear[1:8, ] )

  data.frame( sex=rep( c("male","female","all"), each=16 ),
              transfer=rep( rep(c("ERR","EAR"), each=8), 3 ),
              age_at_exposure=rep(5+0:7*10,3),
              CER=c( rbind( err[,c(1,4,7)], ear[,c(1,4,7)] ) ),
              LowerCI=c( rbind( err[,c(2,5,8)], ear[,c(2,5,8)] ) ),
              UpperCI=c( rbind( err[,c(3,6,9)], err[,c(3,6,9)] ) ),
              CER_all=rep( c( rbind(obj0$err[9,c(1,4,7)],obj0$ear[9,c(1,4,7)]) ), each=8 ),
              LowerCI_all=rep( c( rbind(obj0$err[9,c(2,5,8)],obj0$ear[9,c(2,5,8)]) ), each=8 ),
              UpperCI_all=rep( c( rbind(obj0$err[9,c(3,6,9)],obj0$ear[9,c(3,6,9)]) ), each=8 )
  )
}


plot_CER1 <- function( tab, title="", ylim=c(0,5), drmodel="L" ){  # tab=Mortab$stomach$L_100
  restab1 <- rbind( cbind( get_longtab(tab),  model=drmodel)   )
  restab1$model2 <- paste( restab1$transfer, restab1$model )
  restab1$sex <- factor( c(restab1$sex), levels=c("male","female","all") )
  ggplot( restab1, aes(x=age_at_exposure,y=CER)) + xlab("Age at exposure") + ylab("Cumulatative excess risk (%)") +
    ylim(ylim) + ggtitle(title) +
    facet_wrap(.~sex, scale="free") + geom_line(aes(color=transfer)) +    #linetype=model, )) +
    theme( legend.position=c(.99,.99),
           legend.justification=c("right","top"), legend.box.just="right", legend.margin=margin(6,6,6,6),
           legend.direction = "horizontal", legend.text = element_text(size=8), legend.title = element_text(size=10) ) +
    geom_line(data = restab1, aes(x = age_at_exposure, y = CER_all, color=transfer )) #, linetype=model) )
}

plot_CER2 <- function( tab1, tab2, title="", ylim=c(0,4) ){  # tabL=Mortab$allsolid$L_100; LQtab=Mortab$allsolid$LQ_100
  restab1 <- rbind( cbind( get_longtab(tab1),  model="L" ),
                    cbind( get_longtab(tab2), model="LQ")     )
  restab1$model2 <- paste( restab1$transfer, restab1$model )
  restab1$sex <- factor( c(restab1$sex), levels=c("male","female","all") )
  ggplot( restab1, aes(x=age_at_exposure,y=CER)) + xlab("Age at exposure") + ylab("Cumulatative excess risk (%)") +
    ylim(ylim) + ggtitle(title) +
    facet_wrap(.~sex, scale="free") + geom_line(aes(linetype=model, color=transfer)) +
    theme( legend.position=c(.99,.99),
           legend.justification=c("right","top"), legend.box.just="right", legend.margin=margin(6,6,6,6),
           legend.direction = "horizontal", legend.text = element_text(size=8), legend.title = element_text(size=10) ) +
    geom_line(data = restab1, aes(x = age_at_exposure, y = CER_all, color=transfer, linetype=model) )
}
