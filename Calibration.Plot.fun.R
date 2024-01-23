############################################
############### <calibration ###############
############################################

plot.calibration <- function(PPSdata,alpha=.05,indInt=1,
                             col=customdarkblue,
                             linet=c("solid","dashed"),
                             family="binomial",
                             method = "glm",
                             formula = 'y ~ x',
                             se=F,
                             predx=NULL){
  
  M<-length(PPSdata$PPS)
  
  Models<-c("BB","BB-ED")
  
  data<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    data.frame(y=c(PPSdata$Z[[x]][,1]>qnorm(1-alpha),
                   PPSdata$Z[[x]][,2]>qnorm(1-alpha))*1,
               x=c(PPSdata$PPS[[x]][,indInt+0],
                   PPSdata$PPS[[x]][,indInt+5]))
  }, .parallel = FALSE))
  
  ModelDF<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Models[x],2*nrow(PPSdata$PPS[[x]]))
  }, .parallel = FALSE)
  
  data2=data.frame(data,
                   variable=factor(unlist(ModelDF),
                                   levels=c("BB","BB-ED")))
  
  data2 %>% 
    group_by(variable) %>% 
    do({
      mod = glm(formula = formula,
                family = family,
                data = .)
      mod2 = lm(formula = formula,
                data = .)
      #print(data.frame(Intercept = coef(mod)[1],
      #           Slope = coef(mod)[2]))
      print(round(data.frame(Intercept = coef(mod2)[1],
                             Slope = coef(mod2)[2]),3))
    })
  
  p<-ggplot(data2, aes(x = x, y = y, col=variable, group=variable,linetype=variable)) +
    geom_smooth(method = method,
                method.args = list(family = family),
                formula = formula,
                se=se,
                size=1.25)+
    theme_bw()+labs(colour  = "") +
    scale_color_manual(values=col)+
    #ggtitle(expression(P(Z>z[alpha]~"|"~y[0]~","~y[1])))+
    ggtitle("Calibration curve")+
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 22, 
                                    face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    scale_linetype_manual(name="",values=linet, 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
    theme(legend.position = "none")+
    xlab("predicted probability of positive treatment effect")+ylab("observed average")+
    geom_abline(intercept = 0,slope=1,col="darkgrey",size=1.05)+
    ylim(c(0,1))+xlim(c(0,1))
  
  if(!is.null(predx)){
    pointspred=data.frame(x=predx,
                          y=unlist(plyr::llply(.data = 1:M, .fun = function(x){
                            mod<-glm(formula = formula,
                                     family = family,
                                     data = data2[data2$variable==Models[x],])
                            predict(mod,data.frame(x=predx), type="response")
                          }, .parallel = FALSE)))
    pointspred=pointspred[c(2,4),]
    print(pointspred)
    p <- p + geom_point(pointspred, mapping = 
                          aes(x,y),
                        shape =c(18,16),col=col[c(2,4)],
                        size =7,
                        inherit.aes = FALSE)
  }
  return(p)
}

plot.calibration <- function(PPSuc.1,PPSuc.2,PPSuc.3,alpha=.05,indInt=1,
                             col=c(custompurple,custompurple,customlightblue,customlightblue,customdarkblue,customdarkblue),
                             linet=rep(c("solid","dashed"),3),
                             family="binomial",
                             method = "glm",
                             formula = 'y ~ x',
                             se=F,
                             predx=NULL){
  
  M<-length(PPSuc.1$PPS)
  
  Model1<-c("BB","BB-ED")
  Model2<-c("Logistic","Logistic-ED")
  Model1<-c("NN","NN-ED")
  
  data1<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    data.frame(y=c(PPSuc.1$Z[[x]][,1]>qnorm(1-alpha),
                   PPSuc.1$Z[[x]][,2]>qnorm(1-alpha))*1,
               x=c(PPSuc.1$PPS[[x]][,indInt+0],
                   PPSuc.1$PPS[[x]][,indInt+5]))
  }, .parallel = FALSE))
  data2<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    data.frame(y=c(PPSuc.2$Z[[x]][,1]>qnorm(1-alpha),
                   PPSuc.2$Z[[x]][,2]>qnorm(1-alpha))*1,
               x=c(PPSuc.2$PPS[[x]][,indInt+0],
                   PPSuc.2$PPS[[x]][,indInt+5]))
  }, .parallel = FALSE))
  data3<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    data.frame(y=c(PPSuc.3$Z[[x]][,1]>qnorm(1-alpha),
                   PPSuc.3$Z[[x]][,2]>qnorm(1-alpha))*1,
               x=c(PPSuc.3$PPS[[x]][,indInt+0],
                   PPSuc.3$PPS[[x]][,indInt+5]))
  }, .parallel = FALSE))
  
  ModelDF1<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model1[x],2*nrow(PPSuc.1$PPS[[x]]))
  }, .parallel = FALSE)
  ModelDF2<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model2[x],2*nrow(PPSuc.2$PPS[[x]]))
  }, .parallel = FALSE)
  ModelDF3<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model3[x],2*nrow(PPSuc.3$PPS[[x]]))
  }, .parallel = FALSE)
  
  data_1=data.frame(data1, variable=factor(unlist(ModelDF1), levels=c("BB","BB-ED")))
  data_2=data.frame(data2, variable=factor(unlist(ModelDF2), levels=c("Logistic","Logistic-ED")))
  data_3=data.frame(data3, variable=factor(unlist(ModelDF3), levels=c("NN","NN-ED")))
  
  dataMerge <- rbind(data_1,data_2,data_3)
  
  dataMerge %>% 
    group_by(variable) %>% 
    do({
      mod = glm(formula = formula,
                family = family,
                data = .)
      mod2 = lm(formula = formula,
                data = .)
      #print(data.frame(Intercept = coef(mod)[1],
      #           Slope = coef(mod)[2]))
      print(round(data.frame(Intercept = coef(mod2)[1],
                             Slope = coef(mod2)[2]),3))
    })
  
  p<-ggplot(dataMerge, aes(x = x, y = y, col=variable, group=variable,linetype=variable)) +
    geom_smooth(method = method,
                method.args = list(family = family),
                formula = formula,
                se=se,
                size=1.25)+
    theme_bw()+labs(colour  = "") +
    scale_color_manual(values=col)+
    #ggtitle(expression(P(Z>z[alpha]~"|"~y[0]~","~y[1])))+
    ggtitle("Calibration curve")+
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 22, 
                                    face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    scale_linetype_manual(name="",values=linet, 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
    theme(legend.position = "none")+
    xlab("predicted probability of positive treatment effect")+ylab("observed average")+
    geom_abline(intercept = 0,slope=1,col="darkgrey",size=1.05)+
    ylim(c(0,1))+xlim(c(0,1))
  
  if(!is.null(predx)){
    pointspred=data.frame(x=predx,
                          y=unlist(plyr::llply(.data = 1:M, .fun = function(x){
                            mod<-glm(formula = formula,
                                     family = family,
                                     data = data2[data2$variable==Models[x],])
                            predict(mod,data.frame(x=predx), type="response")
                          }, .parallel = FALSE)))
    pointspred=pointspred[c(2,4),]
    print(pointspred)
    p <- p + geom_point(pointspred, mapping = 
                          aes(x,y),
                        shape =c(18,16),col=col[c(2,4)],
                        size =7,
                        inherit.aes = FALSE)
  }
  return(p)
}

plot.calibration_six <- function(PPSuc.1,PPSuc.2,PPSuc.3,PPSuc.4,PPSuc.5,PPSuc.6,alpha=.05,indInt=1,
                             col=c(black,grey,Rouille,Orange,customlightgreen,customlightblue),
                             linet=rep(c("solid","dotted"),6),
                             family="binomial",
                             method = "glm",
                             formula = 'y ~ x',
                             se=F,
                             predx=NULL){
  
  M<-length(PPSuc.1$PPS)
  
  Model1<-c("BB","BB-ED")
  Model2<-c("Boot","Boot-ED")
  Model3<-c("PB.lo","PB.log-ED")
  Model4<-c("PB.nn","PB.nn-ED")
  Model5<-c("MC.RF","MC.RF-ED")
  Model6<-c("MC.Blog","MC.Blog-ED")
  
  data1<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    data.frame(y=c(PPSuc.1$Z[[x]][,1]>qnorm(1-alpha),
                   PPSuc.1$Z[[x]][,2]>qnorm(1-alpha))*1,
               x=c(PPSuc.1$PPS[[x]][,indInt+0],
                   PPSuc.1$PPS[[x]][,indInt+5]))
  }, .parallel = FALSE))
  data2<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    data.frame(y=c(PPSuc.2$Z[[x]][,1]>qnorm(1-alpha),
                   PPSuc.2$Z[[x]][,2]>qnorm(1-alpha))*1,
               x=c(PPSuc.2$PPS[[x]][,indInt+0],
                   PPSuc.2$PPS[[x]][,indInt+5]))
  }, .parallel = FALSE))
  data3<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    data.frame(y=c(PPSuc.3$Z[[x]][,1]>qnorm(1-alpha),
                   PPSuc.3$Z[[x]][,2]>qnorm(1-alpha))*1,
               x=c(PPSuc.3$PPS[[x]][,indInt+0],
                   PPSuc.3$PPS[[x]][,indInt+5]))
  }, .parallel = FALSE))
  data4<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    data.frame(y=c(PPSuc.4$Z[[x]][,1]>qnorm(1-alpha),
                   PPSuc.4$Z[[x]][,2]>qnorm(1-alpha))*1,
               x=c(PPSuc.4$PPS[[x]][,indInt+0],
                   PPSuc.4$PPS[[x]][,indInt+5]))
  }, .parallel = FALSE))
  data5<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    data.frame(y=c(PPSuc.5$Z[[x]][,1]>qnorm(1-alpha),
                   PPSuc.5$Z[[x]][,2]>qnorm(1-alpha))*1,
               x=c(PPSuc.5$PPS[[x]][,indInt+0],
                   PPSuc.5$PPS[[x]][,indInt+5]))
  }, .parallel = FALSE))
  data6<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    data.frame(y=c(PPSuc.6$Z[[x]][,1]>qnorm(1-alpha),
                   PPSuc.6$Z[[x]][,2]>qnorm(1-alpha))*1,
               x=c(PPSuc.6$PPS[[x]][,indInt+0],
                   PPSuc.6$PPS[[x]][,indInt+5]))
  }, .parallel = FALSE))
  
  
  ModelDF1<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model1[x],2*nrow(PPSuc.1$PPS[[x]]))
  }, .parallel = FALSE)
  ModelDF2<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model2[x],2*nrow(PPSuc.2$PPS[[x]]))
  }, .parallel = FALSE)
  ModelDF3<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model3[x],2*nrow(PPSuc.3$PPS[[x]]))
  }, .parallel = FALSE)
  ModelDF4<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model4[x],2*nrow(PPSuc.4$PPS[[x]]))
  }, .parallel = FALSE)
  ModelDF5<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model5[x],2*nrow(PPSuc.5$PPS[[x]]))
  }, .parallel = FALSE)
  ModelDF6<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model6[x],2*nrow(PPSuc.6$PPS[[x]]))
  }, .parallel = FALSE)
  
  data_1=data.frame(data1, variable=factor(unlist(ModelDF1), levels=c("BB","BB-ED")))
  data_2=data.frame(data2, variable=factor(unlist(ModelDF2), levels=c("Boot","Boot-ED")))
  data_3=data.frame(data3, variable=factor(unlist(ModelDF3), levels=c("PB.log","PB.log-ED")))
  data_4=data.frame(data4, variable=factor(unlist(ModelDF4), levels=c("PB.nn","PB.nn-ED")))
  data_5=data.frame(data5, variable=factor(unlist(ModelDF5), levels=c("MC.RF","MC.RF-ED")))
  data_6=data.frame(data6, variable=factor(unlist(ModelDF6), levels=c("MC.Blog","MC.Blog-ED")))
  
  dataMerge <- rbind(data_1,data_2,data_3,data_4,data_5,data_6)

  dataMerge %>% 
    group_by(variable) %>% 
    do({
      mod = glm(formula = formula,
                family = family,
                data = .)
      mod2 = lm(formula = formula,
                data = .)
      #print(data.frame(Intercept = coef(mod)[1],
      #           Slope = coef(mod)[2]))
      print(round(data.frame(Intercept = coef(mod2)[1],
                             Slope = coef(mod2)[2]),3))
    })
  
  A <- c(rep("BB",800),rep("Boot",800),rep("PB.log",800),rep("PB.nn",800),rep("MC.RF",800),rep("MC.Blog",800))
  line <- rep(c(rep("ID",400),rep("ID-ED",400)),6)
  
  dataMerge$line <- line 
  
  p<-ggplot(dataMerge, aes(x = x, y = y, col=A, group=variable,linetype=line)) +
    geom_smooth(method = method,
                method.args = list(family = family),
                formula = formula,
                se=se,
                size=1.25)+
    theme_bw()+labs(colour  = "black") +
    scale_color_manual(name = "Algorithm",values=c(customlightblue,black,grey,Rouille,Orange,customlightgreen))+
    scale_linetype_manual(name="Data",values=linet,guide = guide_legend(override.aes=aes(fill=NA)))+
    #ggtitle(expression(P(Z>z[alpha]~"|"~y[0]~","~y[1])))+
    ggtitle("Calibration curve")+
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 22, 
                                    face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    theme()+
    xlab("predicted probability of positive treatment effect")+ylab("observed average")+
    geom_abline(intercept = 0,slope=1,col="darkgrey",size=1.05)+
    #ylim(c(0,1))+xlim(c(0,1))
  
  if(!is.null(predx)){
    pointspred=data.frame(x=predx,
                          y=unlist(plyr::llply(.data = 1:M, .fun = function(x){
                            mod<-glm(formula = formula,
                                     family = family,
                                     data = data2[data2$variable==Models[x],])
                            predict(mod,data.frame(x=predx), type="response")
                          }, .parallel = FALSE)))
    pointspred=pointspred[c(2,4),]
    print(pointspred)
    p <- p + geom_point(pointspred, mapping = 
                          aes(x,y),
                        shape =c(18,16),col=col[c(2,4)],
                        size =7,
                        inherit.aes = FALSE)
  }
  return(p)
}
