############### ROC Curve ###############
ROCfun<-function(pH0,pH1,cut=seq(0,1,length=1000)){
  x=apply(t(cut),2,function(x)mean(pH0>x,na.rm = T))
  y=apply(t(cut),2,function(x)mean(pH1>x,na.rm = T))
  out=data.frame(fpr=x,
                 tnr=y)
  return(out)
}

plot.ROC<-function(PPSdata,indInt=1,smoothind=T,alpha=.05,
                   col=c(black,black),
                   linet=c("solid","dashed"),
                   pointROC=NULL){
  
  M<-length(PPSdata$PPS)
  
  Models<-c("BB","BB-ED")
  
  dataROC<-plyr::llply(.data = 1:M, .fun = function(x){
    ROCfun(c(PPSdata$PPS[[x]][,indInt+5][PPSdata$Z[[x]][,2]<qnorm(1-alpha)],
             PPSdata$PPS[[x]][,indInt+0][PPSdata$Z[[x]][,1]<qnorm(1-alpha)]),
           c(PPSdata$PPS[[x]][,indInt+5][PPSdata$Z[[x]][,2]>qnorm(1-alpha)],
             PPSdata$PPS[[x]][,indInt+0][PPSdata$Z[[x]][,1]>qnorm(1-alpha)]))
  }, .parallel = FALSE)
  
  data<-do.call("rbind",dataROC)
  
  pred<-plyr::llply(.data = 1:M, .fun = function(x){
    c(PPSdata$PPS[[x]][,indInt+5][PPSdata$Z[[x]][,2]<qnorm(1-alpha)],
      PPSdata$PPS[[x]][,indInt+0][PPSdata$Z[[x]][,1]<qnorm(1-alpha)],
      PPSdata$PPS[[x]][,indInt+5][PPSdata$Z[[x]][,2]>qnorm(1-alpha)],
      PPSdata$PPS[[x]][,indInt+0][PPSdata$Z[[x]][,1]>qnorm(1-alpha)])
  }, .parallel = FALSE)
  
  obs<-plyr::llply(.data = 1:M, .fun = function(x){
    c(rep(0,length(c(PPSdata$PPS[[x]][,indInt+5][PPSdata$Z[[x]][,2]<qnorm(1-alpha)],
                     PPSdata$PPS[[x]][,indInt+0][PPSdata$Z[[x]][,1]<qnorm(1-alpha)]))),
      rep(1,length(c(PPSdata$PPS[[x]][,indInt+5][PPSdata$Z[[x]][,2]>qnorm(1-alpha)],
                     PPSdata$PPS[[x]][,indInt+0][PPSdata$Z[[x]][,1]>qnorm(1-alpha)]))))
  }, .parallel = FALSE)
  
  dataAUC <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    PRROC::roc.curve(scores.class0 = pred[[x]], weights.class0=obs[[x]],
                     curve=TRUE)$auc
  }, .parallel = FALSE))
  
  # dataAUC<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
  #   DescTools::AUC(dataROC[[x]]$fpr,dataROC[[x]]$tnr)
  # }, .parallel = FALSE))
  
  print(paste("AUC",round(dataAUC,3)))
  
  ModelDF<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Models[x],1000)
  }, .parallel = FALSE)
  
  data2=data.frame(data,
                   variable=factor(unlist(ModelDF),
                                   levels=c("BB","BB-ED")))
  
  p<-ggplot(data2, 
            aes(x=fpr, y=tnr, col=variable, group=variable,linetype=variable)) + 
    geom_abline(intercept = 0,slope=1,col="darkgrey",size=1.05)+
    theme_bw()+labs(colour  = "") +
    scale_color_manual(values=col)+
    #ggtitle(expression(P(Z>z[alpha]~"|"~y[0]~","~y[1])))+
    ggtitle("ROC curve")+
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 22, 
                                    face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    scale_linetype_manual(name="",values=linet, 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
    theme(legend.position = "none")+
    xlab("1-specificity")+ylab("sensitivity")
  
  if(smoothind==T){
    p<-p+geom_smooth(method = "loess", formula = y ~ x,span=.45,
                     size=1.25, se = FALSE)
  }else{
    p<-p+geom_line(size=1.25) 
  }
  
  if(!is.null(pointROC)){
    out=dataROC[[4]][c(1:1000)[seq(0,1,length=1000)>=pointROC][1],]
    print(out)
    p <- p + geom_point(out, mapping = 
                          aes(x=fpr, y=tnr),
                        shape =16,col=col[4],
                        size =7,
                        inherit.aes = FALSE)
    
  }
  
  return(p)
}

plot.ROC.triple <- function(PPSuc.1,PPSuc.2,PPSuc.3,indInt=1,smoothind=T,alpha=.05,
                            col=c(customlightblue,customlightblue,customdarkblue,customdarkblue,custompurple,custompurple),
                            linet=rep(c("solid","dashed"),3),
                            pointROC=NULL){
  M<-length(PPSuc.1$PPS)
  
  Model1<-c("BB","BB-ED")
  Model2<-c("Logistic","Logistic-ED")
  Model3<-c("NN","NN-ED")
  
  dataROC_1<-plyr::llply(.data = 1:M, .fun = function(x){
    ROCfun(c(PPSuc.1$PPS[[x]][,indInt+5][PPSuc.1$Z[[x]][,2]<qnorm(1-alpha)],
             PPSuc.1$PPS[[x]][,indInt+0][PPSuc.1$Z[[x]][,1]<qnorm(1-alpha)]),
           c(PPSuc.1$PPS[[x]][,indInt+5][PPSuc.1$Z[[x]][,2]>qnorm(1-alpha)],
             PPSuc.1$PPS[[x]][,indInt+0][PPSuc.1$Z[[x]][,1]>qnorm(1-alpha)]))
  }, .parallel = FALSE)
  dataROC_2<-plyr::llply(.data = 1:M, .fun = function(x){
    ROCfun(c(PPSuc.2$PPS[[x]][,indInt+5][PPSuc.2$Z[[x]][,2]<qnorm(1-alpha)],
             PPSuc.2$PPS[[x]][,indInt+0][PPSuc.2$Z[[x]][,1]<qnorm(1-alpha)]),
           c(PPSuc.2$PPS[[x]][,indInt+5][PPSuc.2$Z[[x]][,2]>qnorm(1-alpha)],
             PPSuc.2$PPS[[x]][,indInt+0][PPSuc.2$Z[[x]][,1]>qnorm(1-alpha)]))
  }, .parallel = FALSE)
  dataROC_3<-plyr::llply(.data = 1:M, .fun = function(x){
    ROCfun(c(PPSuc.3$PPS[[x]][,indInt+5][PPSuc.3$Z[[x]][,2]<qnorm(1-alpha)],
             PPSuc.3$PPS[[x]][,indInt+0][PPSuc.3$Z[[x]][,1]<qnorm(1-alpha)]),
           c(PPSuc.3$PPS[[x]][,indInt+5][PPSuc.3$Z[[x]][,2]>qnorm(1-alpha)],
             PPSuc.3$PPS[[x]][,indInt+0][PPSuc.3$Z[[x]][,1]>qnorm(1-alpha)]))
  }, .parallel = FALSE)
  
  data1<-do.call("rbind",dataROC_1)
  data2<-do.call("rbind",dataROC_2)
  data3<-do.call("rbind",dataROC_3)
  
  pred1<-plyr::llply(.data = 1:M, .fun = function(x){
    c(PPSuc.1$PPS[[x]][,indInt+5][PPSuc.1$Z[[x]][,2]<qnorm(1-alpha)],
      PPSuc.1$PPS[[x]][,indInt+0][PPSuc.1$Z[[x]][,1]<qnorm(1-alpha)],
      PPSuc.1$PPS[[x]][,indInt+5][PPSuc.1$Z[[x]][,2]>qnorm(1-alpha)],
      PPSuc.1$PPS[[x]][,indInt+0][PPSuc.1$Z[[x]][,1]>qnorm(1-alpha)])
  }, .parallel = FALSE)
  pred2<-plyr::llply(.data = 1:M, .fun = function(x){
    c(PPSuc.2$PPS[[x]][,indInt+5][PPSuc.3$Z[[x]][,2]<qnorm(1-alpha)],
      PPSuc.2$PPS[[x]][,indInt+0][PPSuc.3$Z[[x]][,1]<qnorm(1-alpha)],
      PPSuc.2$PPS[[x]][,indInt+5][PPSuc.3$Z[[x]][,2]>qnorm(1-alpha)],
      PPSuc.2$PPS[[x]][,indInt+0][PPSuc.3$Z[[x]][,1]>qnorm(1-alpha)])
  }, .parallel = FALSE)
  pred3<-plyr::llply(.data = 1:M, .fun = function(x){
    c(PPSuc.3$PPS[[x]][,indInt+5][PPSuc.3$Z[[x]][,2]<qnorm(1-alpha)],
      PPSuc.3$PPS[[x]][,indInt+0][PPSuc.3$Z[[x]][,1]<qnorm(1-alpha)],
      PPSuc.3$PPS[[x]][,indInt+5][PPSuc.3$Z[[x]][,2]>qnorm(1-alpha)],
      PPSuc.3$PPS[[x]][,indInt+0][PPSuc.3$Z[[x]][,1]>qnorm(1-alpha)])
  }, .parallel = FALSE)
   
  obs1<-plyr::llply(.data = 1:M, .fun = function(x){
    c(rep(0,length(c(PPSuc.1$PPS[[x]][,indInt+5][PPSuc.1$Z[[x]][,2]<qnorm(1-alpha)],
                     PPSuc.1$PPS[[x]][,indInt+0][PPSuc.1$Z[[x]][,1]<qnorm(1-alpha)]))),
      rep(1,length(c(PPSuc.1$PPS[[x]][,indInt+5][PPSuc.1$Z[[x]][,2]>qnorm(1-alpha)],
                     PPSuc.1$PPS[[x]][,indInt+0][PPSuc.1$Z[[x]][,1]>qnorm(1-alpha)]))))
  }, .parallel = FALSE)
  obs2<-plyr::llply(.data = 1:M, .fun = function(x){
    c(rep(0,length(c(PPSuc.2$PPS[[x]][,indInt+5][PPSuc.2$Z[[x]][,2]<qnorm(1-alpha)],
                     PPSuc.2$PPS[[x]][,indInt+0][PPSuc.2$Z[[x]][,1]<qnorm(1-alpha)]))),
      rep(1,length(c(PPSuc.2$PPS[[x]][,indInt+5][PPSuc.2$Z[[x]][,2]>qnorm(1-alpha)],
                     PPSuc.2$PPS[[x]][,indInt+0][PPSuc.2$Z[[x]][,1]>qnorm(1-alpha)]))))
  }, .parallel = FALSE)
  obs3<-plyr::llply(.data = 1:M, .fun = function(x){
    c(rep(0,length(c(PPSuc.3$PPS[[x]][,indInt+5][PPSuc.3$Z[[x]][,2]<qnorm(1-alpha)],
                     PPSuc.3$PPS[[x]][,indInt+0][PPSuc.3$Z[[x]][,1]<qnorm(1-alpha)]))),
      rep(1,length(c(PPSuc.3$PPS[[x]][,indInt+5][PPSuc.3$Z[[x]][,2]>qnorm(1-alpha)],
                     PPSuc.3$PPS[[x]][,indInt+0][PPSuc.3$Z[[x]][,1]>qnorm(1-alpha)]))))
  }, .parallel = FALSE)
  
  dataAUC_1 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    PRROC::roc.curve(scores.class0 = pred1[[x]], weights.class0=obs1[[x]],
                     curve=TRUE)$auc
  }, .parallel = FALSE))
  dataAUC_2 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    PRROC::roc.curve(scores.class0 = pred2[[x]], weights.class0=obs2[[x]],
                     curve=TRUE)$auc
  }, .parallel = FALSE))
  dataAUC_3 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    PRROC::roc.curve(scores.class0 = pred3[[x]], weights.class0=obs3[[x]],
                     curve=TRUE)$auc
  }, .parallel = FALSE))
  
  #print(round(dataAUC_1,3))
  
  ModelDF1<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model1[x],1000)
  }, .parallel = FALSE)
  ModelDF2<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model2[x],1000)
  }, .parallel = FALSE)
  ModelDF3<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model3[x],1000)
  }, .parallel = FALSE)
  
  data_1=data.frame(data1,variable=factor(unlist(ModelDF1),levels=c("BB","BB-ED")))
  data_2=data.frame(data2,variable=factor(unlist(ModelDF2),levels=c("Logistic","Logistic-ED")))
  data_3=data.frame(data3,variable=factor(unlist(ModelDF3),levels=c("NN","NN-ED")))
  
  dataMerge <- rbind(data_1,data_2,data_3)
  
  p<-ggplot(dataMerge, 
            aes(x=fpr, y=tnr, col=variable, group=variable,linetype=variable)) + 
    geom_abline(intercept = 0,slope=1,col="darkgrey",size=1.05)+
    theme_bw()+labs(colour  = "") +
    scale_color_manual(values=col)+
    #ggtitle(expression(P(Z>z[alpha]~"|"~y[0]~","~y[1])))+
    ggtitle("ROC curve")+
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 22, 
                                    face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    scale_linetype_manual(name="",values=linet, 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
    theme(legend.position = "none")+
    xlab("1-specificity")+ylab("sensitivity")
  
  if(smoothind==T){
    p<-p+geom_smooth(method = "loess", formula = y ~ x,span=.45,
                     size=1.25, se = FALSE)
  }else{
    p<-p+geom_line(size=1.25) 
  }
  
  if(!is.null(pointROC)){
    out=dataROC[[4]][c(1:1000)[seq(0,1,length=1000)>=pointROC][1],]
    print(out)
    p <- p + geom_point(out, mapping = 
                          aes(x=fpr, y=tnr),
                        shape =16,col=col[4],
                        size =7,
                        inherit.aes = FALSE)
    
  }
  
  return(p)
}

plot.ROC_six <- function(PPSuc.1,PPSuc.2,PPSuc.3,PPSuc.4,PPSuc.5,PPSuc.6,indInt=1,smoothind=T,alpha=.05,
                            col=c(customlightblue,black,grey,Rouille,Orange,customlightgreen),
                            linet=rep(c("solid","dotted"),6),
                            pointROC=NULL){
  M<-length(PPSuc.1$PPS)
  
  Model1<-c("BB","BB-ED")
  Model2<-c("Boot","Boot-ED")
  Model3<-c("PB.log","PB.log-ED")
  Model4<-c("PB.nn","PB.nn-ED")
  Model5<-c("MC.RF","MC.RF-ED")
  Model6<-c("MC.Blog","MC.Blog-ED")
  
  dataROC_1<-plyr::llply(.data = 1:M, .fun = function(x){
    ROCfun(c(PPSuc.1$PPS[[x]][,indInt+5][PPSuc.1$Z[[x]][,2]<qnorm(1-alpha)],
             PPSuc.1$PPS[[x]][,indInt+0][PPSuc.1$Z[[x]][,1]<qnorm(1-alpha)]),
           c(PPSuc.1$PPS[[x]][,indInt+5][PPSuc.1$Z[[x]][,2]>qnorm(1-alpha)],
             PPSuc.1$PPS[[x]][,indInt+0][PPSuc.1$Z[[x]][,1]>qnorm(1-alpha)]))
  }, .parallel = FALSE)
  dataROC_2<-plyr::llply(.data = 1:M, .fun = function(x){
    ROCfun(c(PPSuc.2$PPS[[x]][,indInt+5][PPSuc.2$Z[[x]][,2]<qnorm(1-alpha)],
             PPSuc.2$PPS[[x]][,indInt+0][PPSuc.2$Z[[x]][,1]<qnorm(1-alpha)]),
           c(PPSuc.2$PPS[[x]][,indInt+5][PPSuc.2$Z[[x]][,2]>qnorm(1-alpha)],
             PPSuc.2$PPS[[x]][,indInt+0][PPSuc.2$Z[[x]][,1]>qnorm(1-alpha)]))
  }, .parallel = FALSE)
  dataROC_3<-plyr::llply(.data = 1:M, .fun = function(x){
    ROCfun(c(PPSuc.3$PPS[[x]][,indInt+5][PPSuc.3$Z[[x]][,2]<qnorm(1-alpha)],
             PPSuc.3$PPS[[x]][,indInt+0][PPSuc.3$Z[[x]][,1]<qnorm(1-alpha)]),
           c(PPSuc.3$PPS[[x]][,indInt+5][PPSuc.3$Z[[x]][,2]>qnorm(1-alpha)],
             PPSuc.3$PPS[[x]][,indInt+0][PPSuc.3$Z[[x]][,1]>qnorm(1-alpha)]))
  }, .parallel = FALSE)
  dataROC_4<-plyr::llply(.data = 1:M, .fun = function(x){
    ROCfun(c(PPSuc.4$PPS[[x]][,indInt+5][PPSuc.4$Z[[x]][,2]<qnorm(1-alpha)],
             PPSuc.4$PPS[[x]][,indInt+0][PPSuc.4$Z[[x]][,1]<qnorm(1-alpha)]),
           c(PPSuc.4$PPS[[x]][,indInt+5][PPSuc.4$Z[[x]][,2]>qnorm(1-alpha)],
             PPSuc.4$PPS[[x]][,indInt+0][PPSuc.4$Z[[x]][,1]>qnorm(1-alpha)]))
  }, .parallel = FALSE)
  dataROC_5<-plyr::llply(.data = 1:M, .fun = function(x){
    ROCfun(c(PPSuc.5$PPS[[x]][,indInt+5][PPSuc.5$Z[[x]][,2]<qnorm(1-alpha)],
             PPSuc.5$PPS[[x]][,indInt+0][PPSuc.5$Z[[x]][,1]<qnorm(1-alpha)]),
           c(PPSuc.5$PPS[[x]][,indInt+5][PPSuc.5$Z[[x]][,2]>qnorm(1-alpha)],
             PPSuc.5$PPS[[x]][,indInt+0][PPSuc.5$Z[[x]][,1]>qnorm(1-alpha)]))
  }, .parallel = FALSE)
  dataROC_6<-plyr::llply(.data = 1:M, .fun = function(x){
    ROCfun(c(PPSuc.6$PPS[[x]][,indInt+5][PPSuc.6$Z[[x]][,2]<qnorm(1-alpha)],
             PPSuc.6$PPS[[x]][,indInt+0][PPSuc.6$Z[[x]][,1]<qnorm(1-alpha)]),
           c(PPSuc.6$PPS[[x]][,indInt+5][PPSuc.6$Z[[x]][,2]>qnorm(1-alpha)],
             PPSuc.6$PPS[[x]][,indInt+0][PPSuc.6$Z[[x]][,1]>qnorm(1-alpha)]))
  }, .parallel = FALSE)
  
  data1<-do.call("rbind",dataROC_1)
  data2<-do.call("rbind",dataROC_2)
  data3<-do.call("rbind",dataROC_3)
  data4<-do.call("rbind",dataROC_4)
  data5<-do.call("rbind",dataROC_5)
  data6<-do.call("rbind",dataROC_6)
  
  pred1<-plyr::llply(.data = 1:M, .fun = function(x){
    c(PPSuc.1$PPS[[x]][,indInt+5][PPSuc.1$Z[[x]][,2]<qnorm(1-alpha)],
      PPSuc.1$PPS[[x]][,indInt+0][PPSuc.1$Z[[x]][,1]<qnorm(1-alpha)],
      PPSuc.1$PPS[[x]][,indInt+5][PPSuc.1$Z[[x]][,2]>qnorm(1-alpha)],
      PPSuc.1$PPS[[x]][,indInt+0][PPSuc.1$Z[[x]][,1]>qnorm(1-alpha)])
  }, .parallel = FALSE)
  pred2<-plyr::llply(.data = 1:M, .fun = function(x){
    c(PPSuc.2$PPS[[x]][,indInt+5][PPSuc.2$Z[[x]][,2]<qnorm(1-alpha)],
      PPSuc.2$PPS[[x]][,indInt+0][PPSuc.2$Z[[x]][,1]<qnorm(1-alpha)],
      PPSuc.2$PPS[[x]][,indInt+5][PPSuc.2$Z[[x]][,2]>qnorm(1-alpha)],
      PPSuc.2$PPS[[x]][,indInt+0][PPSuc.2$Z[[x]][,1]>qnorm(1-alpha)])
  }, .parallel = FALSE)
  pred3<-plyr::llply(.data = 1:M, .fun = function(x){
    c(PPSuc.3$PPS[[x]][,indInt+5][PPSuc.3$Z[[x]][,2]<qnorm(1-alpha)],
      PPSuc.3$PPS[[x]][,indInt+0][PPSuc.3$Z[[x]][,1]<qnorm(1-alpha)],
      PPSuc.3$PPS[[x]][,indInt+5][PPSuc.3$Z[[x]][,2]>qnorm(1-alpha)],
      PPSuc.3$PPS[[x]][,indInt+0][PPSuc.3$Z[[x]][,1]>qnorm(1-alpha)])
  }, .parallel = FALSE)
  pred4<-plyr::llply(.data = 1:M, .fun = function(x){
    c(PPSuc.4$PPS[[x]][,indInt+5][PPSuc.4$Z[[x]][,2]<qnorm(1-alpha)],
      PPSuc.4$PPS[[x]][,indInt+0][PPSuc.4$Z[[x]][,1]<qnorm(1-alpha)],
      PPSuc.4$PPS[[x]][,indInt+5][PPSuc.4$Z[[x]][,2]>qnorm(1-alpha)],
      PPSuc.4$PPS[[x]][,indInt+0][PPSuc.4$Z[[x]][,1]>qnorm(1-alpha)])
  }, .parallel = FALSE)
  pred5<-plyr::llply(.data = 1:M, .fun = function(x){
    c(PPSuc.5$PPS[[x]][,indInt+5][PPSuc.5$Z[[x]][,2]<qnorm(1-alpha)],
      PPSuc.5$PPS[[x]][,indInt+0][PPSuc.5$Z[[x]][,1]<qnorm(1-alpha)],
      PPSuc.5$PPS[[x]][,indInt+5][PPSuc.5$Z[[x]][,2]>qnorm(1-alpha)],
      PPSuc.5$PPS[[x]][,indInt+0][PPSuc.5$Z[[x]][,1]>qnorm(1-alpha)])
  }, .parallel = FALSE)
  pred6<-plyr::llply(.data = 1:M, .fun = function(x){
    c(PPSuc.6$PPS[[x]][,indInt+5][PPSuc.6$Z[[x]][,2]<qnorm(1-alpha)],
      PPSuc.6$PPS[[x]][,indInt+0][PPSuc.6$Z[[x]][,1]<qnorm(1-alpha)],
      PPSuc.6$PPS[[x]][,indInt+5][PPSuc.6$Z[[x]][,2]>qnorm(1-alpha)],
      PPSuc.6$PPS[[x]][,indInt+0][PPSuc.6$Z[[x]][,1]>qnorm(1-alpha)])
  }, .parallel = FALSE)
  
  obs1<-plyr::llply(.data = 1:M, .fun = function(x){
    c(rep(0,length(c(PPSuc.1$PPS[[x]][,indInt+5][PPSuc.1$Z[[x]][,2]<qnorm(1-alpha)],
                     PPSuc.1$PPS[[x]][,indInt+0][PPSuc.1$Z[[x]][,1]<qnorm(1-alpha)]))),
      rep(1,length(c(PPSuc.1$PPS[[x]][,indInt+5][PPSuc.1$Z[[x]][,2]>qnorm(1-alpha)],
                     PPSuc.1$PPS[[x]][,indInt+0][PPSuc.1$Z[[x]][,1]>qnorm(1-alpha)]))))
  }, .parallel = FALSE)
  obs2<-plyr::llply(.data = 1:M, .fun = function(x){
    c(rep(0,length(c(PPSuc.2$PPS[[x]][,indInt+5][PPSuc.2$Z[[x]][,2]<qnorm(1-alpha)],
                     PPSuc.2$PPS[[x]][,indInt+0][PPSuc.2$Z[[x]][,1]<qnorm(1-alpha)]))),
      rep(1,length(c(PPSuc.2$PPS[[x]][,indInt+5][PPSuc.2$Z[[x]][,2]>qnorm(1-alpha)],
                     PPSuc.2$PPS[[x]][,indInt+0][PPSuc.2$Z[[x]][,1]>qnorm(1-alpha)]))))
  }, .parallel = FALSE)
  obs3<-plyr::llply(.data = 1:M, .fun = function(x){
    c(rep(0,length(c(PPSuc.3$PPS[[x]][,indInt+5][PPSuc.3$Z[[x]][,2]<qnorm(1-alpha)],
                     PPSuc.3$PPS[[x]][,indInt+0][PPSuc.3$Z[[x]][,1]<qnorm(1-alpha)]))),
      rep(1,length(c(PPSuc.3$PPS[[x]][,indInt+5][PPSuc.3$Z[[x]][,2]>qnorm(1-alpha)],
                     PPSuc.3$PPS[[x]][,indInt+0][PPSuc.3$Z[[x]][,1]>qnorm(1-alpha)]))))
  }, .parallel = FALSE)
  obs4<-plyr::llply(.data = 1:M, .fun = function(x){
    c(rep(0,length(c(PPSuc.4$PPS[[x]][,indInt+5][PPSuc.4$Z[[x]][,2]<qnorm(1-alpha)],
                     PPSuc.4$PPS[[x]][,indInt+0][PPSuc.4$Z[[x]][,1]<qnorm(1-alpha)]))),
      rep(1,length(c(PPSuc.4$PPS[[x]][,indInt+5][PPSuc.4$Z[[x]][,2]>qnorm(1-alpha)],
                     PPSuc.4$PPS[[x]][,indInt+0][PPSuc.4$Z[[x]][,1]>qnorm(1-alpha)]))))
  }, .parallel = FALSE)
  obs5<-plyr::llply(.data = 1:M, .fun = function(x){
    c(rep(0,length(c(PPSuc.5$PPS[[x]][,indInt+5][PPSuc.5$Z[[x]][,2]<qnorm(1-alpha)],
                     PPSuc.5$PPS[[x]][,indInt+0][PPSuc.5$Z[[x]][,1]<qnorm(1-alpha)]))),
      rep(1,length(c(PPSuc.5$PPS[[x]][,indInt+5][PPSuc.5$Z[[x]][,2]>qnorm(1-alpha)],
                     PPSuc.5$PPS[[x]][,indInt+0][PPSuc.5$Z[[x]][,1]>qnorm(1-alpha)]))))
  }, .parallel = FALSE)
  obs6<-plyr::llply(.data = 1:M, .fun = function(x){
    c(rep(0,length(c(PPSuc.6$PPS[[x]][,indInt+5][PPSuc.6$Z[[x]][,2]<qnorm(1-alpha)],
                     PPSuc.6$PPS[[x]][,indInt+0][PPSuc.6$Z[[x]][,1]<qnorm(1-alpha)]))),
      rep(1,length(c(PPSuc.6$PPS[[x]][,indInt+5][PPSuc.6$Z[[x]][,2]>qnorm(1-alpha)],
                     PPSuc.6$PPS[[x]][,indInt+0][PPSuc.6$Z[[x]][,1]>qnorm(1-alpha)]))))
  }, .parallel = FALSE)
  
  dataAUC_1 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    PRROC::roc.curve(scores.class0 = pred1[[x]], weights.class0=obs1[[x]],
                     curve=TRUE)$auc
  }, .parallel = FALSE))
  dataAUC_2 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    PRROC::roc.curve(scores.class0 = pred2[[x]], weights.class0=obs2[[x]],
                     curve=TRUE)$auc
  }, .parallel = FALSE))
  dataAUC_3 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    PRROC::roc.curve(scores.class0 = pred3[[x]], weights.class0=obs3[[x]],
                     curve=TRUE)$auc
  }, .parallel = FALSE))
  dataAUC_4 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    PRROC::roc.curve(scores.class0 = pred4[[x]], weights.class0=obs4[[x]],
                     curve=TRUE)$auc
  }, .parallel = FALSE))
  dataAUC_5 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    PRROC::roc.curve(scores.class0 = pred5[[x]], weights.class0=obs5[[x]],
                     curve=TRUE)$auc
  }, .parallel = FALSE))
  dataAUC_6 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    PRROC::roc.curve(scores.class0 = pred6[[x]], weights.class0=obs6[[x]],
                     curve=TRUE)$auc
  }, .parallel = FALSE))
  
  #print(round(dataAUC_1,3))
  
  ModelDF1<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model1[x],1000)
  }, .parallel = FALSE)
  ModelDF2<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model2[x],1000)
  }, .parallel = FALSE)
  ModelDF3<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model3[x],1000)
  }, .parallel = FALSE)
  ModelDF4<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model4[x],1000)
  }, .parallel = FALSE)
  ModelDF5<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model5[x],1000)
  }, .parallel = FALSE)
  ModelDF6<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model6[x],1000)
  }, .parallel = FALSE)
  
  data_1=data.frame(data1,variable=factor(unlist(ModelDF1),levels=c("BB","BB-ED")))
  data_2=data.frame(data2,variable=factor(unlist(ModelDF2),levels=c("Boot","Boot-ED")))
  data_3=data.frame(data3,variable=factor(unlist(ModelDF3),levels=c("PB.log","PB.log-ED")))
  data_4=data.frame(data4,variable=factor(unlist(ModelDF4),levels=c("PB.nn","PB.nn-ED")))
  data_5=data.frame(data5,variable=factor(unlist(ModelDF5),levels=c("MC.RF","MC.RF-ED")))
  data_6=data.frame(data6,variable=factor(unlist(ModelDF6),levels=c("MC.Blog","MC.Blog-ED")))
  
  dataMerge <- rbind(data_1,data_2,data_3,data_4,data_5,data_6)
  A <- c(rep("BB",2000),rep("Boot",2000),rep("PB.log",2000),rep("PB.nn",2000),rep("MC.RF",2000),rep("MC.Blog",2000))
  line <- rep(c(rep("ID",1000),rep("ID-ED",1000)),6)
  
  p<-ggplot(dataMerge, 
            aes(x=fpr, y=tnr, col=A, group=variable,linetype=line)) + 
    geom_abline(intercept = 0,slope=1,col="darkgrey",size=1.05)+
    theme_bw()+labs(colour  = "") +
    scale_color_manual(name = "Algorithm",values=c(customlightblue,black,grey,Rouille,Orange,customlightgreen))+
    scale_linetype_manual(name="Data",values=linet,guide = guide_legend(override.aes=aes(fill=NA)))+
    #ggtitle(expression(P(Z>z[alpha]~"|"~y[0]~","~y[1])))+
    ggtitle("ROC curve")+
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 22, 
                                    face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    theme()+
    xlab("1-specificity")+ylab("sensitivity")
  
#  if(smoothind==T){
#    p<-p+geom_smooth(method = "loess", formula = y ~ x,span=.45,
#                     size=1.25, se = FALSE)
#  }else{
    p<-p+geom_line(size=1.25) 
#  }
  
  if(!is.null(pointROC)){
    out=dataROC[[4]][c(1:1000)[seq(0,1,length=1000)>=pointROC][1],]
    print(out)
    p <- p + geom_point(out, mapping = 
                          aes(x=fpr, y=tnr),
                        shape =16,col=col[4],
                        size =7,
                        inherit.aes = FALSE)
    
  }
  
  return(p)
}
