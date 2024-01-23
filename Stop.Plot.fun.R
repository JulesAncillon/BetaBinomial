####################################
############### Stop ###############
####################################

plot.stop <- function(PPSdata,alpha=.05,indInt=1,smoothind=T,
                      cut=seq(0,1,length=1000),futility=T,
                      col=c(customgreen,customgreen),
                      linet=c("solid","dashed"),
                      pointb=NULL){
  
  M<-length(PPSdata$PPS)
  
  Models<-c("BB","BB-ED")
  
  data<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    if(futility==T){
      y=apply(t(cut),2,function(k)mean(PPSdata$PPS[[x]][,indInt+5]<k,na.rm = T))
    }else{
      y=apply(t(cut),2,function(k)mean(PPSdata$PPS[[x]][,indInt+0]>k,na.rm = T))
    }
    data.frame(y=y,
               x=cut)
  }, .parallel = FALSE))
  
  ModelDF<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Models[x],length(cut))
  }, .parallel = FALSE)
  
  data2=data.frame(data,
                   variable=factor(unlist(ModelDF),
                                   levels=c("BB","BB-ED")))
  
  p<-ggplot(data2, aes(x = x, y = y, col=variable, group=variable,linetype=variable)) +
    theme_bw()+labs(colour  = "") +
    scale_color_manual(values=col)+
    ggtitle("Probability of futility stopping when there is no treatment effect")+
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 19.5, 
                                    face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    scale_linetype_manual(name="",values=linet, 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
    theme_bw()+labs(colour  = "black") +
    xlab("bound")+ylab("probability")+
    ylim(c(0,1))+xlim(c(0,1))
  
  if(smoothind==T){
    p<-p+geom_smooth(method = "loess", formula = y ~ x,span=.075,
                     size=1.25, se = FALSE)
  }else{
    p<-p+geom_line(size=1.25) 
  }
  
  if(!is.null(pointb)){
    dataM=data2[data2$variable=="LR-ED",]
    out = dataM[c(1:1000)[seq(0,1,length=1000)>=pointb][1],]
    print(out)
    p <- p + geom_point(out, mapping = 
                          aes(x, y),
                        shape =16,col=col[4],
                        size =7,
                        inherit.aes = FALSE)
  }
  return(p)
}


plot.triple.stop <- function(PPSuc.1,PPSuc.2,PPSuc.3,alpha=.05,indInt=1,smoothind=T,
                      cut=seq(0,1,length=1000),futility=T,
                      col=c('#999999','#999999',"#1D2951","#1D2951","#BFAD87","#BFAD87"),
                      linet=rep(c("solid","dashed"),3),
                      pointb=NULL){
  
  M<-length(PPSuc.1$PPS)
  
  Model1<-c("BB","BB-ED")
  Model2<-c("PB_LOG","PB_LOG-ED")
  Model3<-c("PB_NN","PB_NN-ED")
  
  data1 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    if(futility==T){
      y=apply(t(cut),2,function(k)mean(PPSuc.1$PPS[[x]][,indInt+5]<k,na.rm = T))
    }else{
      y=apply(t(cut),2,function(k)mean(PPSuc.1$PPS[[x]][,indInt+0]>k,na.rm = T))
    }
    data.frame(y=y,
               x=cut)
  }, .parallel = FALSE))
  data2 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    if(futility==T){
      y=apply(t(cut),2,function(k)mean(PPSuc.2$PPS[[x]][,indInt+5]<k,na.rm = T))
    }else{
      y=apply(t(cut),2,function(k)mean(PPSuc.2$PPS[[x]][,indInt+0]>k,na.rm = T))
    }
    data.frame(y=y,
               x=cut)
  }, .parallel = FALSE))
  data3 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    if(futility==T){
      y=apply(t(cut),2,function(k)mean(PPSuc.3$PPS[[x]][,indInt+5]<k,na.rm = T))
    }else{
      y=apply(t(cut),2,function(k)mean(PPSuc.3$PPS[[x]][,indInt+0]>k,na.rm = T))
    }
    data.frame(y=y,
               x=cut)
  }, .parallel = FALSE))
  
  ModelDF1<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model1[x],length(cut))
  }, .parallel = FALSE)
  ModelDF2<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model2[x],length(cut))
  }, .parallel = FALSE)
  ModelDF3<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model3[x],length(cut))
  }, .parallel = FALSE)
  
  data_1 =data.frame(data1, variable=factor(unlist(ModelDF1),levels=c("BB","BB-ED")))
  data_2 =data.frame(data2, variable=factor(unlist(ModelDF2),levels=c("PB_LOG","PB_LOG-ED")))
  data_3 =data.frame(data3, variable=factor(unlist(ModelDF3),levels=c("PB_NN","PB_NN-ED")))
  
  dataMerge <- rbind(data_1,data_2,data_3)
  
  p<-ggplot(dataMerge, aes(x = x, y = y, col=variable, group=variable,linetype=variable)) +
    theme_bw()+labs(colour  = "") +
    scale_color_manual(values=col)+
    ggtitle("Probability of futility stopping when there is no treatment effect")+
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 19.5, 
                                    face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    scale_linetype_manual(name="",values=linet, 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
    theme(legend.position = "none")+
    xlab("bound")+ylab("probability")+
    ylim(c(0,1))+xlim(c(0,1))
  
  if(smoothind==T){
    p<-p+geom_smooth(method = "loess", formula = y ~ x,span=.075,
                     size=1, se = FALSE)
  }else{
    p<-p+geom_line(size=1) 
  }
  
  if(!is.null(pointb)){
    dataM=dataMerge[dataMerge$variable=="LR-ED",]
    out = dataM[c(1:1000)[seq(0,1,length=1000)>=pointb][1],]
    print(out)
    p <- p + geom_point(out, mapping = 
                          aes(x, y),
                        shape =16,col=col[4],
                        size =6,
                        inherit.aes = FALSE)
  }
  return(p)
}



plot.six.stop <- function(PPSuc.1,PPSuc.2,PPSuc.3,PPSuc.4,PPSuc.5,PPSuc.6,alpha=.05,indInt=1,smoothind=T,
                             cut=seq(0,1,length=1000),futility=T,
                             col=c(black,grey,Rouille,Orange,customlightgreen,customlightblue),
                             linet=rep(c("solid","dotted"),6),
                             pointb=NULL){
  
  M<-length(PPSuc.1$PPS)
  
  Model1<-c("BB","BB-ED")
  Model2<-c("Boot","Boot-ED")
  Model3<-c("PB.log","PB.log-ED")
  Model4<-c("PB.nn","PB.nn-ED")
  Model5<-c("MC.RF","MC.RF-ED")
  Model6<-c("MC.Blog","MC.Blog-ED")
  
  data1 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    if(futility==T){
      y=apply(t(cut),2,function(k)mean(PPSuc.1$PPS[[x]][,indInt+5]<k,na.rm = T))
    }else{
      y=apply(t(cut),2,function(k)mean(PPSuc.1$PPS[[x]][,indInt+0]>k,na.rm = T))
    }
    data.frame(y=y,
               x=cut)
  }, .parallel = FALSE))
  data2 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    if(futility==T){
      y=apply(t(cut),2,function(k)mean(PPSuc.2$PPS[[x]][,indInt+5]<k,na.rm = T))
    }else{
      y=apply(t(cut),2,function(k)mean(PPSuc.2$PPS[[x]][,indInt+0]>k,na.rm = T))
    }
    data.frame(y=y,
               x=cut)
  }, .parallel = FALSE))
  data3 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    if(futility==T){
      y=apply(t(cut),2,function(k)mean(PPSuc.3$PPS[[x]][,indInt+5]<k,na.rm = T))
    }else{
      y=apply(t(cut),2,function(k)mean(PPSuc.3$PPS[[x]][,indInt+0]>k,na.rm = T))
    }
    data.frame(y=y,
               x=cut)
  }, .parallel = FALSE))
  data4 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    if(futility==T){
      y=apply(t(cut),2,function(k)mean(PPSuc.4$PPS[[x]][,indInt+5]<k,na.rm = T))
    }else{
      y=apply(t(cut),2,function(k)mean(PPSuc.4$PPS[[x]][,indInt+0]>k,na.rm = T))
    }
    data.frame(y=y,
               x=cut)
  }, .parallel = FALSE))
  data5 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    if(futility==T){
      y=apply(t(cut),2,function(k)mean(PPSuc.5$PPS[[x]][,indInt+5]<k,na.rm = T))
    }else{
      y=apply(t(cut),2,function(k)mean(PPSuc.5$PPS[[x]][,indInt+0]>k,na.rm = T))
    }
    data.frame(y=y,
               x=cut)
  }, .parallel = FALSE))
  data6 <-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    if(futility==T){
      y=apply(t(cut),2,function(k)mean(PPSuc.6$PPS[[x]][,indInt+5]<k,na.rm = T))
    }else{
      y=apply(t(cut),2,function(k)mean(PPSuc.6$PPS[[x]][,indInt+0]>k,na.rm = T))
    }
    data.frame(y=y,
               x=cut)
  }, .parallel = FALSE))
  
  ModelDF1<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model1[x],length(cut))
  }, .parallel = FALSE)
  ModelDF2<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model2[x],length(cut))
  }, .parallel = FALSE)
  ModelDF3<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model3[x],length(cut))
  }, .parallel = FALSE)
  ModelDF4<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model4[x],length(cut))
  }, .parallel = FALSE)
  ModelDF5<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model5[x],length(cut))
  }, .parallel = FALSE)
  ModelDF6<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Model6[x],length(cut))
  }, .parallel = FALSE)
  
  data_1 =data.frame(data1, variable=factor(unlist(ModelDF1),levels=c("BB","BB-ED")))
  data_2 =data.frame(data2, variable=factor(unlist(ModelDF2),levels=c("Boot","Boot-ED")))
  data_3 =data.frame(data3, variable=factor(unlist(ModelDF3),levels=c("PB.log","PB.log-ED")))
  data_4 =data.frame(data4, variable=factor(unlist(ModelDF4),levels=c("PB.nn","PB.nn-ED")))
  data_5 =data.frame(data5, variable=factor(unlist(ModelDF5),levels=c("MC.RF","MC.RF-ED")))
  data_6 =data.frame(data6, variable=factor(unlist(ModelDF6),levels=c("MC.Blog","MC.Blog-ED")))
  
  A <- c(rep("BB",2000),rep("Boot",2000),rep("PB.log",2000),rep("PB.nn",2000),rep("MC.RF",2000),rep("MC.Blog",2000))
  line <- rep(c(rep("ID",1000),rep("ID-ED",1000)),6)
  
  dataMerge <- rbind(data_1,data_2,data_3,data_4,data_5,data_6)
  
  p<-ggplot(dataMerge, aes(x = x, y = y, col=A, group=variable,linetype=line)) +
    theme_bw()+labs(colour  = "") +
    scale_color_manual(name = "Algorithm",values=c(customlightblue,black,grey,Rouille,Orange,customlightgreen))+
    ggtitle("Probability of true futility stopping")+
    scale_linetype_manual(name="Data",values=linet, 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 22, 
                                    face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    xlab("bound")+ylab("Probability of futility stopping when there is no treatment effect")+
    ylim(c(0,1))+xlim(c(0,1))
  
  if(smoothind==T){
    p<-p+geom_smooth(method = "loess", formula = y ~ x,span=.075,
                     size=1, se = FALSE)
  }else{
    p<-p+geom_line(size=1) 
  }
  
  if(!is.null(pointb)){
    dataM=dataMerge[dataMerge$variable=="LR-ED",]
    out = dataM[c(1:1000)[seq(0,1,length=1000)>=pointb][1],]
    print(out)
    p <- p + geom_point(out, mapping = 
                          aes(x, y),
                        shape =16,col=col[4],
                        size =6,
                        inherit.aes = FALSE)
  }
  return(p)
}
