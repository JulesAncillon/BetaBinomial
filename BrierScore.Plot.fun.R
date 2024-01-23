library(ggplot2)

############### Brier Score ###############
BSf<-function(PPSd,Z,alpha=.05){
  c((PPSd-(Z>qnorm(1-alpha)*1))^2)
}

data.BS<-function(PPS,alpha=.05,n.int,interims=5){
  M<-length(PPS$PPS)
  data<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPS$PPS[[x]][,k],PPS$Z[[x]][,1],alpha),
             BSf(PPS$PPS[[x]][,k+5],PPS$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  rownames(data) <- c("BB","BB-ED")
  colnames(data) <- n.int[1:5]
  return(data)
}


plot.BS<-function(dataBS,ylim=c(0,.5), col=c(Black,Black),
                  linet=c("solid","dashed")){ 
  dataBS2<-reshape2::melt(dataBS)
  
  p<- ggplot(dataBS2, 
             aes(x=Var2, y=value, col=Var1, group=Var1,linetype=Var1)) + 
    geom_point(shape=16,size=3)+
    geom_line(size=1.25) +
    theme_bw()+labs(colour  = "") +
    scale_color_manual(values=col)+
    #ggtitle(expression(P(Z>z[alpha]~"|"~y[0]~","~y[1])))+
    #ylim(ylim)+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    xlab("number of patients enrolled before interim")+ylab("Brier score")+
    scale_linetype_manual(name="",values=linet, 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
    theme(legend.position = "none")
  show(p)
  return(dataBS2)
}

######################################
# Double
######################################

data.BS_double<-function(PPS_1, PPS_2,alpha=.05,n.int,interims=5){
  M<-length(PPS_1$PPS)
  data_1<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPS_1$PPS[[x]][,k],PPS_1$Z[[x]][,1],alpha),
             BSf(PPS_1$PPS[[x]][,k+5],PPS_1$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_2<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPS_2$PPS[[x]][,k],PPS_2$Z[[x]][,1],alpha),
             BSf(PPS_2$PPS[[x]][,k+5],PPS_2$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  rownames(data_1) <- c("1L-10N-ID","1L-10N-ED")
  colnames(data_1) <- n.int
  rownames(data_2) <- c("10L-10N-ID","10L-10N-ED")
  colnames(data_2) <- n.int
  return(list(data_1,data_2))
}

plot.BS_double<-function(dataBS_double,ylim=c(0,.35),
                         linet=rep(c("solid","dashed"),3)){ # Solid = ID & dashed = ED
  
  dataBS_1<-reshape2::melt(dataBS_double[[1]])
  dataBS_2<-reshape2::melt(dataBS_double[[2]])
  
  dataBS <- rbind(dataBS_1,dataBS_2)
  
  p <- ggplot(dataBS, 
              aes(x=Var2, y=value, col=Var1, group=Var1,linetype=Var1)) + 
    geom_point(shape=16,size=3)+ geom_line(size=1.25) +theme_bw() + labs(colour  = "") +
    scale_color_manual(values=c('#885717','#885717','#03C04A' ,'#03C04A' ))+
    #ylim(ylim)+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    xlab("number of patients enrolled before interim")+ylab("Brier score")+
    scale_linetype_manual(name="",values=linet, 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
    theme(legend.position = "top")
  show(p)
  #return(dataBS2)
}


######################################
# Triple
######################################

data.BS_triple<-function(PPSuc.1, PPSuc.2, PPSuc.3, alpha=.05,n.int,interims=5){
  M<-length(PPSuc.1$PPS)
 
  data_1<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.1$PPS[[x]][,k],PPSuc.1$Z[[x]][,1],alpha),
             BSf(PPSuc.1$PPS[[x]][,k+5],PPSuc.1$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_2<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.2$PPS[[x]][,k],PPSuc.2$Z[[x]][,1],alpha),
             BSf(PPSuc.2$PPS[[x]][,k+5],PPSuc.2$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_3<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.3$PPS[[x]][,k],PPSuc.3$Z[[x]][,1],alpha),
             BSf(PPSuc.3$PPS[[x]][,k+5],PPSuc.3$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  rownames(data_1) <- c("BB-ID","BB-ID+ED")
  colnames(data_1) <- n.int
  rownames(data_2) <- c("PB.GLM-ID","PB.GLM-ID+ED")
  colnames(data_2) <- n.int
  rownames(data_3) <- c("PB.NN-ID","PB.NN-ID+ED")
  colnames(data_3) <- n.int
  return(list(data_1,data_2,data_3))
}

plot.BS_triple <- function(dataBS_triple, ylim = c(0, .45), linet = rep(c("solid", "dotted"), 4)) {
  dataBS_1 <- reshape2::melt(dataBS_triple[[1]])
  dataBS_2 <- reshape2::melt(dataBS_triple[[2]])
  dataBS_3 <- reshape2::melt(dataBS_triple[[3]])
  
  dataBS <- rbind(dataBS_1, dataBS_2, dataBS_3)
  line <- rep(c("ID", "ID+ED"), 15)
  
  ggplot(dataBS,
         aes(x = Var2, y = value, col = Var1, group = Var1, linetype = line)) +
    geom_point(aes(shape = Var1), size = 3) +  # Added shape aesthetics here
    geom_line(size = 1.25) +
    theme_bw() +
    labs(colour = "") +
    scale_color_manual(name = "Algorithm", values = c('#b2df8a', '#b2df8a', '#a6cee3', '#a6cee3', '#1f78b4', '#1f78b4')) +
    ylim(.0, .35) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size = 30)) +
    xlab("number of patients enrolled before interim") + ylab("Brier score") +
    scale_linetype_manual(name = "", values = linet,
                          guide = guide_legend(override.aes = aes(fill = NA))) +
    scale_shape_manual(name = "Algorithm", values = c(1,1,2, 2,3, 3)) +  # Define shapes here
    theme(legend.position = "top")
}



######################################
# Four
######################################
data.BS_four<-function(PPSuc.0, PPSuc.1, PPSuc.2, PPSuc.3, alpha=.05,n.int,interims=5){
  M<-length(PPSuc.0$PPS)
  data_0<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.0$PPS[[x]][,k],PPSuc.0$Z[[x]][,1],alpha),
             BSf(PPSuc.0$PPS[[x]][,k+5],PPSuc.0$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_1<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.1$PPS[[x]][,k],PPSuc.1$Z[[x]][,1],alpha),
             BSf(PPSuc.1$PPS[[x]][,k+5],PPSuc.1$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_2<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.2$PPS[[x]][,k],PPSuc.2$Z[[x]][,1],alpha),
             BSf(PPSuc.2$PPS[[x]][,k+5],PPSuc.2$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_3<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.3$PPS[[x]][,k],PPSuc.3$Z[[x]][,1],alpha),
             BSf(PPSuc.3$PPS[[x]][,k+5],PPSuc.3$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  rownames(data_0) <- c("Classic","Classic-ED")
  colnames(data_0) <- n.int
  rownames(data_1) <- c("Rand.Ratio","Rand.Ratio-ED")
  colnames(data_1) <- n.int
  rownames(data_2) <- c("Biomarkers","Biomarkers-ED")
  colnames(data_2) <- n.int
  rownames(data_3) <- c("Outcome","Outcome-ED")
  colnames(data_3) <- n.int
  return(list(data_0,data_1,data_2,data_3))
}


plot.BS_four<-function(dataBS_four,ylim=c(0,.20),
                         linet=rep(c("solid","dashed"),4)){ # Solid = ID & dashed = ED
  
  dataBS_0<-reshape2::melt(dataBS_four[[1]])
  dataBS_1<-reshape2::melt(dataBS_four[[2]])
  dataBS_2<-reshape2::melt(dataBS_four[[3]])
  dataBS_3<-reshape2::melt(dataBS_four[[4]])
  
  dataBS <- rbind(dataBS_0,dataBS_1,dataBS_2,dataBS_3)
  
  ggplot(dataBS, 
         aes(x=Var2, y=value, col=Var1, group=Var1,linetype=Var1)) + 
    geom_point(shape=16,size=3)+ geom_line(size=1.25) +theme_bw() + labs(colour  = "") +
    scale_color_manual(values=c(Black, Black, customlightblue,customlightblue, custompurple,custompurple,customdarkblue,customdarkblue))+
    ylim(ylim)+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    xlab("number of patients enrolled before interim")+ylab("Brier score")+
    scale_linetype_manual(name="",values=linet, 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
    theme(legend.position = "top")
}


######################################
# Five
######################################

data.BS_five<-function(PPSuc.0, PPSuc.1, PPSuc.2, PPSuc.3, PPSuc.4, alpha=.05,n.int,interims=5){
  
  if ((length(PPSuc.0$PPS) != length(PPSuc.1$PPS)) | (length(PPSuc.0$PPS) != length(PPSuc.2$PPS)) | (length(PPSuc.0$PPS) != length(PPSuc.3$PPS)) | (length(PPSuc.0$PPS) != length(PPSuc.4$PPS))) {
    return("Erro Length")
  }
  M<-length(PPSuc.0$PPS)
  
  data_0<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.0$PPS[[x]][,k],PPSuc.0$Z[[x]][,1],alpha),
             BSf(PPSuc.0$PPS[[x]][,k+5],PPSuc.0$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_1<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.1$PPS[[x]][,k],PPSuc.1$Z[[x]][,1],alpha),
             BSf(PPSuc.1$PPS[[x]][,k+5],PPSuc.1$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_2<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.2$PPS[[x]][,k],PPSuc.2$Z[[x]][,1],alpha),
             BSf(PPSuc.2$PPS[[x]][,k+5],PPSuc.2$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_3<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.3$PPS[[x]][,k],PPSuc.3$Z[[x]][,1],alpha),
             BSf(PPSuc.3$PPS[[x]][,k+5],PPSuc.3$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_4<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.4$PPS[[x]][,k],PPSuc.4$Z[[x]][,1],alpha),
             BSf(PPSuc.4$PPS[[x]][,k+5],PPSuc.4$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  
  rownames(data_0) <- c("BB","BB-ED")
  colnames(data_0) <- n.int
  rownames(data_1) <- c("Boot","Boot-ED")
  colnames(data_1) <- n.int
  rownames(data_2) <- c("PB.log","PB.log-ED")
  colnames(data_2) <- n.int
  rownames(data_3) <- c("PB.nn","PB.nn-ED")
  colnames(data_3) <- n.int
  rownames(data_3) <- c("PB.nn","PB.nn-ED")
  colnames(data_3) <- n.int
  rownames(data_4) <- c("MC.RF","MC.RF-ED")
  colnames(data_4) <- n.int
  return(list(data_0,data_1,data_2,data_3,data_4))
}

plot.BS_five<-function(dataBS_six,ylim=c(0,.35),
                      linet=rep(c("solid","dotted"),5)){ # Solid = ID & dashed = ED
  
  dataBS_0<-reshape2::melt(dataBS_six[[1]])
  dataBS_1<-reshape2::melt(dataBS_six[[2]])
  dataBS_2<-reshape2::melt(dataBS_six[[3]])
  dataBS_3<-reshape2::melt(dataBS_six[[4]])
  dataBS_4<-reshape2::melt(dataBS_six[[5]])
  
  dataBS <- rbind(dataBS_0,dataBS_1,dataBS_2,dataBS_3,dataBS_4)
  A <- c(rep("BB",10),rep("Boot",10),rep("PB.log",10),rep("PB.nn",10),rep("MC.RF",10))
  col <- c(rep(Black,10), rep(grey,10), rep(customlightblue,10), rep(Rouille,10),rep(Orange,10))
  
  ggplot(dataBS, 
         aes(x=Var2, y=value, col= A, group=Var1,linetype=rep(c("ID","ID-ED"),30))) + 
    geom_point(shape=16,size=3)+ geom_line(size=1.25) +theme_bw() + labs(colour  = "") +
    scale_color_manual(name = "Algorithm",values=c(customlightblue,Black,grey,Rouille,Orange,customlightgreen))+
    scale_linetype_manual(name="Data",values=linet,guide = guide_legend(override.aes=aes(fill=NA)))+
    #ylim(ylim)+
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 22, 
                                    face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    xlab("number of patients enrolled before interim")+ylab("Brier score") +
    ggtitle("Brier Score")+
    labs(colour  = "Black")+
    theme()
}

######################################
# Six
######################################
data.BS_six<-function(PPSuc.0, PPSuc.1, PPSuc.2, PPSuc.3, PPSuc.4, PPSuc.5, alpha=.05,n.int,interims=5){
  if ((length(PPSuc.0$PPS) != length(PPSuc.1$PPS)) | (length(PPSuc.0$PPS) != length(PPSuc.2$PPS)) | (length(PPSuc.0$PPS) != length(PPSuc.3$PPS)) | (length(PPSuc.0$PPS) != length(PPSuc.4$PPS))| (length(PPSuc.0$PPS) != length(PPSuc.5$PPS))) {
    return("Erro Length")
  }
  M<-length(PPSuc.0$PPS)
  
  data_0<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.0$PPS[[x]][,k],PPSuc.0$Z[[x]][,1],alpha),
             BSf(PPSuc.0$PPS[[x]][,k+5],PPSuc.0$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_1<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.1$PPS[[x]][,k],PPSuc.1$Z[[x]][,1],alpha),
             BSf(PPSuc.1$PPS[[x]][,k+5],PPSuc.1$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_2<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.2$PPS[[x]][,k],PPSuc.2$Z[[x]][,1],alpha),
             BSf(PPSuc.2$PPS[[x]][,k+5],PPSuc.2$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_3<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.3$PPS[[x]][,k],PPSuc.3$Z[[x]][,1],alpha),
             BSf(PPSuc.3$PPS[[x]][,k+5],PPSuc.3$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_4<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.4$PPS[[x]][,k],PPSuc.4$Z[[x]][,1],alpha),
             BSf(PPSuc.4$PPS[[x]][,k+5],PPSuc.4$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  data_5<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    sapply(1:5,function(k) 
      mean(c(BSf(PPSuc.5$PPS[[x]][,k],PPSuc.5$Z[[x]][,1],alpha),
             BSf(PPSuc.5$PPS[[x]][,k+5],PPSuc.5$Z[[x]][,2],alpha))))
    
  }, .parallel = FALSE))
  
  rownames(data_0) <- c("Reg_00","Reg_00-ED")
  colnames(data_0) <- n.int
  rownames(data_1) <- c("Reg_02","Reg_02-ED")
  colnames(data_1) <- n.int
  rownames(data_2) <- c("Reg_04","Reg_04-ED")
  colnames(data_2) <- n.int
  rownames(data_3) <- c("Reg_06","Reg_06-ED")
  colnames(data_3) <- n.int
  rownames(data_4) <- c("Reg_08","Reg_08-ED")
  colnames(data_4) <- n.int
  rownames(data_5) <- c("Reg_10","Reg_10-ED")
  colnames(data_5) <- n.int
  return(list(data_0,data_1,data_2,data_3,data_4,data_5))
}

plot.BS_six<-function(dataBS_six,ylim=c(0,.35),
                       linet=rep(c("solid","dotted"),6)){ # Solid = ID & dashed = ED
  
  dataBS_0<-reshape2::melt(dataBS_six[[1]])
  dataBS_1<-reshape2::melt(dataBS_six[[2]])
  dataBS_2<-reshape2::melt(dataBS_six[[3]])
  dataBS_3<-reshape2::melt(dataBS_six[[4]])
  dataBS_4<-reshape2::melt(dataBS_six[[5]])
  dataBS_5<-reshape2::melt(dataBS_six[[6]])
  
  dataBS <- rbind(dataBS_0,dataBS_1,dataBS_2,dataBS_3,dataBS_4,dataBS_5)
  A <- c(rep("Reg_00",10),rep("Reg_02",10),rep("Reg_04",10),rep("Reg_06",10),rep("Reg_08",10),rep("Reg_10",10))
  col <- c(rep("black",10), rep("79828D",10), rep("blue",10), rep("985717",10),rep("orange",10),rep("90EE90",10))
  
  ggplot(dataBS, 
         aes(x=Var2, y=value, col= A, group=Var1,linetype=rep(c("ID","ID-ED"),30))) + 
    geom_point(shape=16,size=3)+ geom_line(size=1.25) +theme_bw() + labs(colour  = "") +
    scale_color_manual(name = "Algorithm",values=c("blue","black","79828D","985717","orange","90EE90"))+
    scale_linetype_manual(name="Data",values=linet,guide = guide_legend(override.aes=aes(fill=NA)))+
    #ylim(ylim)+
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 22, 
                                    face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    xlab("number of patients enrolled before interim")+ylab("Brier score") +
    ggtitle("Brier Score")+
    labs(colour  = "Black")+
    theme()
}



plot.BSS<-function(BS.NN, BS.BB,ylim=c(0,.35),
                   col=c(customgreen,customgreen),
                   linet=rep(c("solid","dashed"),3)){ # Solid = ID & dashed = ED
  
  BS.BB_ID <- BS.BB[BS.BB$Var1 == "BB",]
  BS.BB_ED <- BS.BB[BS.BB$Var1 != "BB",]
  
  BS.NN_ID <- BS.NN[BS.NN$Var1 == "BB",]
  BS.NN_ED <- BS.NN[BS.NN$Var1 != "BB",]
  
  data.BS_ID <- BS.NN_ID
  data.BS_ID[,3] <- 1 - BS.NN_ID[,3]/BS.BB_ID[,3]
  data.BS_ED = BS.NN_ED 
  data.BS_ED[,3] <- 1 - BS.NN_ED[,3]/BS.BB_ED[,3]
  
  plot(x= data.BS_ED[,2], y = c(-1,rep(1,length(data.BS_ED[,3])-1)), type = "l", col = "white", xlab = "Number of patients enrolled before interim", ylab = "% of improvement",
       main=paste("Brier Skill Score"))
  points(x= data.BS_ID[,2], y = data.BS_ID[,3], type = "l", col = "Black")
  points(x= data.BS_ED[,2], y = data.BS_ED[,3], type = "b", col = "Black")
  return(data.frame(data.BS_ID[,3],data.BS_ED[,3]))
}


