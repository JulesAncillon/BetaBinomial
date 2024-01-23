plot.true.stop <- function(PPSdata,bound=.1,futility=T,alpha=.05,
                            n.int,
                            col=c('#999999','#999999'),
                            linet=c("solid","dashed"),
                            family="binomial",
                            method = "glm",
                            formula = 'y ~ x',
                            se=F){
  
  M<-length(PPSdata$PPS)
  
  data<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    if(futility==T){
      sapply(1:5,function(k) 
        mean(PPSdata$PPS[[x]][,k+5][PPSdata$Z[[x]][,1]>qnorm(1-alpha)]<bound))
    }else{
      sapply(1:5,function(k) 
        mean(PPSdata$PPS[[x]][,k+5][PPSdata$Z[[x]][,2]<qnorm(1-alpha)]>bound))
    }
  }, .parallel = FALSE))
  
  rownames(data) <- c("BB","BB-ED")
  colnames(data) <- n.int
  
  print(round(data,3))
  
  data2<-reshape2::melt(data)
  
  ggplot(data2, 
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
    xlab("number of patients enrolled before interim")+ylab("discontinuation probability (positive treatment)")+
    scale_linetype_manual(name="",values=linet, 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
    theme(legend.position = "none")
  
}













plot.true.stop_triple <- function(PPSuc.1, PPSuc.2, PPSuc.3, bound = .1, futility = T, alpha = .05,
                                  n.int,
                                  col = c('#b2df8a', '#a6cee3', '#1f78b4'),
                                  linet = rep(c("solid", "dotted"), 3),
                                  family = "binomial",
                                  method = "glm",
                                  formula = 'y ~ x',
                                  se = F) {
  
  M <- length(PPSuc.1$PPS)
  
  Model1 <- c("BB-ID", "BB-ID+ED")
  Model2 <- c("PB_GLM-ID", "PB_LOG-ID+ED")
  Model3 <- c("PB_NN-ID", "PB_NN-ID+ED")
  
  data1 <- do.call("rbind", plyr::llply(.data = 1:M, .fun = function(x) {
    if (futility == T) {
      sapply(1:5, function(k) 
        mean(PPSuc.1$PPS[[x]][, k + 5][PPSuc.1$Z[[x]][, 1] > qnorm(1 - alpha)] < bound))
    } else {
      sapply(1:5, function(k) 
        mean(PPSuc.1$PPS[[x]][, k + 5][PPSuc.1$Z[[x]][, 2] < qnorm(1 - alpha)] > bound))
    }
  }, .parallel = FALSE))
  data2 <- do.call("rbind", plyr::llply(.data = 1:M, .fun = function(x) {
    if (futility == T) {
      sapply(1:5, function(k) 
        mean(PPSuc.2$PPS[[x]][, k + 5][PPSuc.2$Z[[x]][, 1] > qnorm(1 - alpha)] < bound))
    } else {
      sapply(1:5, function(k) 
        mean(PPSuc.2$PPS[[x]][, k + 5][PPSuc.2$Z[[x]][, 2] < qnorm(1 - alpha)] > bound))
    }
  }, .parallel = FALSE))
  data3 <- do.call("rbind", plyr::llply(.data = 1:M, .fun = function(x) {
    if (futility == T) {
      sapply(1:5, function(k) 
        mean(PPSuc.3$PPS[[x]][, k + 5][PPSuc.3$Z[[x]][, 1] > qnorm(1 - alpha)] < bound))
    } else {
      sapply(1:5, function(k) 
        mean(PPSuc.3$PPS[[x]][, k + 5][PPSuc.3$Z[[x]][, 2] < qnorm(1 - alpha)] > bound))
    }
  }, .parallel = FALSE))
  
  rownames(data1) <- c("BB-ID", "BB-ID+ED")
  colnames(data1) <- n.int
  rownames(data2) <- c("PB_GLM-ID", "PB_GLM-ID+ED")
  colnames(data2) <- n.int
  rownames(data3) <- c("PB_NN-ID", "PB_NN-ID+ED")
  colnames(data3) <- n.int
  
  data_1 <- reshape2::melt(data1)
  data_2 <- reshape2::melt(data2)
  data_3 <- reshape2::melt(data3)
  
  dataMerge <- rbind(data_1, data_2, data_3)
  A <- c(rep("BB", 10), rep("PB.GLM", 10), rep("PB.NN", 10))
  line <- rep(c("ID", "ID+ED"), 15)  # Le 30 ici est peu comprehensible
  
  ggplot(dataMerge, 
         aes(x = Var2, y = value, col = A, group = Var1, linetype = line, shape = A)) +
    geom_point(size = 3) +  # Added shape aesthetics here
    geom_line(size = 1.25) +
    theme_bw() +
    labs(colour = "") +
    ggtitle("Probability of True Stop") +
    scale_color_manual(name = "Algorithm", values = col) +
    scale_shape_manual(name = "Algorithm", values = c(1, 2, 3)) +  # Define shapes here
    ylim(c(0, 1)) +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 22,
                                    face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size = 30)) +
    xlab("number of patients enrolled before interim") + ylab("discontinuation probability (positive treatment)") +
    scale_linetype_manual(name = "Data", values = linet,
                          guide = guide_legend(override.aes = aes(fill = NA))) +
    theme(legend.position = "bottom")
}


