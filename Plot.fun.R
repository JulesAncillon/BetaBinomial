#########PLOTS FUN ########
data.PPS<-function(ZH0,ZH1,TrialsH0,TrialsH1,n.int){
  PPS<-list()
  Z<-list()
  
  auxH1 <- do.call("cbind",plyr::llply(.data = 1:length(TrialsH1), .fun = function(x){
    TrialsH1[[x]]$PPS.ID
  }, .parallel = FALSE))
  auxH0 <- do.call("cbind",plyr::llply(.data = 1:length(TrialsH0), .fun = function(x){
    TrialsH0[[x]]$PPS.ID
  }, .parallel = FALSE))
  PPS[[1]] <- cbind(auxH1,auxH0)
  Z[[1]] <- cbind(ZH1,ZH0)
  
  auxH1 <- do.call("cbind",plyr::llply(.data = 1:length(TrialsH1), .fun = function(x){
    TrialsH1[[x]]$PPS.ED
  }, .parallel = FALSE))
  auxH0 <- do.call("cbind",plyr::llply(.data = 1:length(TrialsH0), .fun = function(x){
    TrialsH0[[x]]$PPS.ED
  }, .parallel = FALSE))
  PPS[[2]] <- cbind(auxH1,auxH0)
  Z[[2]] <- Z[[1]]
  
  return(list(PPS=PPS,
              Z=Z))
}

plot.HistLT<-function(PPSdata,panels=2,indInt=1,alpha=.05,titlelist=NULL,
                      col=c(customgreen,customgreen),
                      colH=rep("white",2),
                      linet=c("solid","dashed"),
                      interims=5,patients=166){
  
  interimaux<-round(patients/(interims+1)*c(1:interims))
  
  if(is.null(titlelist)){
    titlelist<-c(paste("No positive treatment effect detected, H0 not rejected"),#, m =",interimaux[indInt]),
                 paste("Positive treatment effect detected, H0 not rejected"),#, m =",interimaux[indInt]),
                 paste("No positive treatment effect detected, H0 rejected"),#, m =",interimaux[indInt]),
                 paste("Positive treatment effect detected, H0 rejected"))#, m =",interimaux[indInt]))
  }
  
  H1<-HistLT(PPSdata,indInt,F,T,titlelist[1],alpha,col,colH,linet)
  H2<-HistLT(PPSdata,indInt,T,T,titlelist[2],alpha,col,colH,linet)
  H3<-HistLT(PPSdata,indInt,F,F,titlelist[3],alpha,col,colH,linet)
  H4<-HistLT(PPSdata,indInt,T,F,titlelist[4],alpha,col,colH,linet)
  
  if(panels==2){
    HistLTp<-gridExtra::grid.arrange(H1, H4,  ncol=2, nrow =1,
                                     top = textGrob("Frequency plots",gp=gpar(fontsize=22,font=2)))
  }else{
    HistLTp<-gridExtra::grid.arrange(H1, H2, H3, H4, ncol=2, nrow =2,
                                     top = textGrob("Frequency plots",gp=gpar(fontsize=22,font=2)))
  }
  
  return(HistLTp)
}

############### TypeIerror ###############
TypeIerror <- function(TrialsH0,n.int, Z, t=0.5, color = 'blue'){
  x = 1:length(n.int)
  thresh = t
  Error1.ID = plyr::llply(.data = 1:length(x), .fun =  function(k){
    length(TrialsH0[[k]]$PPS.ID[TrialsH0[[k]]$PPS.ID>thresh])/length(TrialsH0[[k]]$PPS.ID)
  })
  Error1.ED = plyr::llply(.data = 1:length(x), .fun =  function(k){
    length(TrialsH0[[k]]$PPS.ED[TrialsH0[[k]]$PPS.ED>thresh])/length(TrialsH0[[k]]$PPS.ED)
  })
  l = list(unlist(Error1.ID),unlist(Error1.ED))
  plot(c(0,rep(1,length(n.int)-1)), x = n.int, col = 'white', type ='l', ylab ="Probability of Type I error", xlab = "Number of patients enrolled before interim",
       main=paste("Proportion of Type I errors for t=",thresh))
  points(y = unlist(Error1.ID), x = n.int,type='l',col= color)
  points(rep(1,length(n.int))*Z,x=n.int,  col = 'grey', type = 'l')
  points(y = unlist(Error1.ED), x = n.int,type='b',col= color)
  return(l)
}
TypeIerror.triple <- function(TrialsH0.bb,TrialsH0.log,TrialsH0.nn,n.int,t=0.5){
  x = 1:length(n.int)
  thresh = t
  Errorbb.ID = plyr::llply(.data = 1:length(x), .fun =  function(k){
    length(TrialsH0.bb[[k]]$PPS.ID[TrialsH0.bb[[k]]$PPS.ID>thresh])/length(TrialsH0.bb[[k]]$PPS.ID)
  })
  Errorbb.ED = plyr::llply(.data = 1:length(x), .fun =  function(k){
    length(TrialsH0.bb[[k]]$PPS.ED[TrialsH0.bb[[k]]$PPS.ED>thresh])/length(TrialsH0.bb[[k]]$PPS.ED)
  })
  Errorlog.ID = plyr::llply(.data = 1:length(x), .fun =  function(k){
    length(TrialsH0.log[[k]]$PPS.ID[TrialsH0.log[[k]]$PPS.ID>thresh])/length(TrialsH0.log[[k]]$PPS.ID)
  })
  Errorlog.ED = plyr::llply(.data = 1:length(x), .fun =  function(k){
    length(TrialsH0.log[[k]]$PPS.ED[TrialsH0.log[[k]]$PPS.ED>thresh])/length(TrialsH0.log[[k]]$PPS.ED)
  })
  Errornn.ID = plyr::llply(.data = 1:length(x), .fun =  function(k){
    length(TrialsH0.nn[[k]]$PPS.ID[TrialsH0.nn[[k]]$PPS.ID>thresh])/length(TrialsH0.nn[[k]]$PPS.ID)
  })
  Errornn.ED = plyr::llply(.data = 1:length(x), .fun =  function(k){
    length(TrialsH0.nn[[k]]$PPS.ED[TrialsH0.nn[[k]]$PPS.ED>thresh])/length(TrialsH0.nn[[k]]$PPS.ED)
  })
  l = list(unlist(Errorbb.ID),unlist(Errorbb.ED))
  plot(c(0,rep(1,length(n.int)-1)), x = n.int, col = 'white', type ='l', ylab ="Probability of Type I error", xlab = "Number of patients enrolled before interim",
       main=paste("Proportion of Type I errors for t=",thresh))
  points(rep(1,length(n.int))*0.05,x=n.int,  col = 'grey', type = 'l')
  points(y = unlist(Errorbb.ID), x = n.int,type='l',col= custompurple,lwd=2.0)
  points(y = unlist(Errorbb.ED), x = n.int,type='b',col= custompurple,lwd=2.0)
  points(y = unlist(Errorlog.ID), x = n.int,type='l',col= customlightblue,lwd=2.0)
  points(y = unlist(Errorlog.ED), x = n.int,type='b',col= customlightblue,lwd=2.0)
  points(y = unlist(Errornn.ID), x = n.int,type='l',col= customdarkblue,lwd=2.0)
  points(y = unlist(Errornn.ED), x = n.int,type='b',col= customdarkblue,lwd=2.0)
}

############### Power ###############
Power <- function(TrialsH1,n.int, Z, t = 0.5){
  x = 1:length(n.int)
  thresh = t
  Power.ID = plyr::llply(.data = 1:length(x), .fun =  function(k){
    1-length(TrialsH1[[k]]$PPS.ID[TrialsH1[[k]]$PPS.ID<thresh])/length(TrialsH1[[k]]$PPS.ID)
  })
  Power.ED = plyr::llply(.data = 1:length(x), .fun =  function(k){
    1-length(TrialsH1[[k]]$PPS.ED[TrialsH1[[k]]$PPS.ED<thresh])/length(TrialsH1[[k]]$PPS.ED)
  })
  l = list(unlist(Power.ID),unlist(Power.ED))
  plot(y = c(0,rep(1,length(n.int)-1)),x = n.int, type='l', col = 'white', ylab ="Power", xlab = "Number of patients enrolled before interim",
       main=paste("Statistical Power for t=",thresh))
  points(y = rep(1,length(n.int))*Z,x = n.int, type='l', col = 'grey')
  points(y = unlist(Power.ID), x = n.int,type='l',col='red')
  points(y = unlist(Power.ED), x = n.int,type='b',col='red')
  return(l)
}
Power.triple <- function(TrialsH0.bb,TrialsH0.log,TrialsH0.nn, n.int, t = 0.5){
  x = 1:length(n.int)
  thresh = t
  Power.bb.ID = plyr::llply(.data = 1:length(x), .fun =  function(k){
    1-length(TrialsH0.bb[[k]]$PPS.ID[TrialsH0.bb[[k]]$PPS.ID<thresh])/length(TrialsH0.bb[[k]]$PPS.ID)
  })
  Power.bb.ED = plyr::llply(.data = 1:length(x), .fun =  function(k){
    1-length(TrialsH0.bb[[k]]$PPS.ED[TrialsH0.bb[[k]]$PPS.ED<thresh])/length(TrialsH0.bb[[k]]$PPS.ED)
  })
  Power.log.ID = plyr::llply(.data = 1:length(x), .fun =  function(k){
    1-length(TrialsH0.log[[k]]$PPS.ID[TrialsH0.log[[k]]$PPS.ID<thresh])/length(TrialsH0.log[[k]]$PPS.ID)
  })
  Power.log.ED = plyr::llply(.data = 1:length(x), .fun =  function(k){
    1-length(TrialsH0.log[[k]]$PPS.ED[TrialsH0.log[[k]]$PPS.ED<thresh])/length(TrialsH0.log[[k]]$PPS.ED)
  })
  Power.nn.ID = plyr::llply(.data = 1:length(x), .fun =  function(k){
    1-length(TrialsH0.nn[[k]]$PPS.ID[TrialsH0.nn[[k]]$PPS.ID<thresh])/length(TrialsH0.nn[[k]]$PPS.ID)
  })
  Power.nn.ED = plyr::llply(.data = 1:length(x), .fun =  function(k){
    1-length(TrialsH0.nn[[k]]$PPS.ED[TrialsH0.nn[[k]]$PPS.ED<thresh])/length(TrialsH0.nn[[k]]$PPS.ED)
  })
  plot(y = c(.2,rep(1,length(n.int)-1)),x = n.int, type='l', col = 'white', ylab ="Power", xlab = "Number of patients enrolled before interim",
       main=paste("Statistical Power for t=",thresh))
  points(y = rep(1,length(n.int))*0.05,x = n.int, type='l', col = 'grey')
  points(y = unlist(Power.bb.ID), x = n.int,type='l',col=custompurple, lwd=2.0)
  points(y = unlist(Power.bb.ED), x = n.int,type='b',col=custompurple, lwd=2.0)
  points(y = unlist(Power.log.ID), x = n.int,type='l',col=customlightblue, lwd=2.0)
  points(y = unlist(Power.log.ED), x = n.int,type='b',col=customlightblue, lwd=2.0)
  points(y = unlist(Power.nn.ID), x = n.int,type='l',col=customdarkblue, lwd=2.0)
  points(y = unlist(Power.nn.ED), x = n.int,type='b',col=customdarkblue, lwd=2.0)
}

