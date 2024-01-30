
#' Compute a standard z-test for covariance
#' @param y possible outcome
#' @param a arm
#' @export
#' @example
#' Outcome = rbinom(100,1,.30)
#' Arm = sample(rep(0:1,100/2),replace=FALSE)
#' z.test(y=Outcome,a=Arm)
z.test<-function(y,a){
  patients<-c(sum(a==0),sum(a!=0))
  outcome<-c(sum(y[a==0]),sum(y[a!=0]))
  p.success<-outcome/patients
  p.t<-sum(outcome)/sum(patients)
  delta<-p.success[2]-p.success[1]
  t.std<-delta/sqrt(p.t*(1-p.t)/patients[1]+p.t*(1-p.t)/patients[2])
  # return(list(z.t.statistic=t.std,
  #             z.te=delta))
  return(t.std)
}

z.test2<-function(XC,XT,NC,NT){
  patients<-c(NC,NT)
  outcome<-c(XC,XT)
  p.success<-outcome/patients
  p.t<-sum(outcome)/sum(patients)
  delta<-p.success[2]-p.success[1]
  t.std<-delta/sqrt(p.t*(1-p.t)/patients[1]+p.t*(1-p.t)/patients[2])
  #t.std<-delta/sqrt(p.t*(1-p.t)*(1/patients[1]+1/patients[2]))
  #return(list(z.ts=t.std,
  #            z.te=delta))
  return(t.std)
}

#' Generate possible outcomes
#' @param N The number of patients in the patiential outcome dataset
#' @param p.X A vector p of probabilities for the covariates
#' @return An N times p matrix of response outcome.
#' @export
#' @examples
#' X.sample(N=100, p.X=c(.25, .5, .9))
X.sample<-function(N, p.X, cont = FALSE, sigma = .1){
  N.x = length(p.X)
  X=1*(matrix(runif(N.x*N), N, N.x) <= matrix(p.X, N, N.x, byrow=T))
  names <- rep(" ",length(p.X))
  
  if (cont == TRUE) {
    for (i in 1:N.x) {
      X[,i] = rnorm(N[1],p.X[i],sigma)
    }
  }
  
  for (i in 1:length(p.X)) { names[i] <- paste("X",i,sep="")}
  colnames(X) <- names
  
  return(X)
}

#' Generate possible outcomes of internal trial
#' @param N A vector with the maximum number of patients in each arm of the patiential outcome dataset
#' @param p.X A vector p of probabilities for the covariates
#' @param lambda Arrival rate of patients
#' @param beta.arm coefficient of arm
#' @param beta.x coefficient of covariates
#' @export
#' @examples
#' sim.simple.trial.arm(N=rep(100,2),p.X=c(.5,.5,.5,.5,.5),lambda=10,beta.arm=c(-0.5,0),beta.x=c(0.14,-0.15,0.65,0.2,0.1))
sim.simple.trial.arm<-function(N,p.X,lambda,beta.arm,beta.x, interaction = F, cont=F){
  N.all   = sum(N)
  k       = length(beta.arm)-1
  X.A     = sample(apply(as.matrix(0:k),2,function(x)rep(x,N[x+1])),replace=FALSE)
  X       = X.sample(sum(N),p.X,cont)
  Y <- output(X,beta.arm,beta.x,interaction)
  Y_hat<-do.call("rbind",plyr::llply(.data = 1:N.all, .fun = function(x){
    Y[x,X.A[x]+1]
  }, .parallel = FALSE))
  colnames(Y) = paste0("Y.A", 0:k)
  arr.time= c(0,cumsum(rexp(n=(N.all-1),  rate=lambda)))
  
  return(data.frame(arr.time=arr.time,
                    X.A=X.A,
                    X,
                    ID=1,
                    Y,
                    Y=Y_hat))
}

# X = X.ED
# beta.arm = c(0,0,qnorm(.7)-qnorm(.5))
# beta.x = c(qnorm(.7),-1.5,1,0)
# interaction = F

output <- function(X,beta.arm,beta.x,interaction=F){
  if (interaction == T) {
    Y       = t(apply(X,1,function(x) runif(length(beta.arm)) <= pnorm(beta.arm + c(beta.x[1] + x%*%beta.x[-1] + 2*x[1]*x[3] + 10*x[1]*x[2]*x[3] )) ))*1
  } 
  else  { 
    Y       = t(apply(X,1,function(x)runif(length(beta.arm))<=pnorm(beta.arm + c(beta.x[1] + x%*%beta.x[-1]))))*1
  }
} 

###################################################################################
###################################- PLOTS FUN -###################################
###################################################################################

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
  colnames(data) <- n.int
  return(data)
}

plot.BS<-function(dataBS,ylim=c(0,.35),
                  col=c(customgreen,customgreen),
                  linet=rep(c("solid","dashed"),3)){
  
  dataBS2<-reshape2::melt(dataBS)
  
  ggplot(dataBS2, 
         aes(x=Var2, y=value, col=Var1, group=Var1,linetype=Var1)) + 
    geom_point(shape=16,size=3)+
    geom_line(size=1.25) +
    theme_bw()+labs(colour  = "") +
    scale_color_manual(values=col)+
    #ggtitle(expression(P(Z>z[alpha]~"|"~y[0]~","~y[1])))+
    ylim(ylim)+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.title = element_text(size =30))+
    xlab("number of patients enrolled before interim")+ylab("Brier score")+
    scale_linetype_manual(name="",values=linet, 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
    theme(legend.position = "none")
}

HistLT<-function(PPSdata,indInt,TP=TRUE,H0=TRUE,titlep,alpha=.05,
                 col=c(customgreen,customgreen),
                 colH=rep("white",2),
                 linet=c("solid","dashed"),
                 bw=0.05){
  M<-length(PPSdata$PPS)
  indH0<-ifelse(H0==T,2,1)
  indPPSH0<-ifelse(H0==T,5,0)
  Models<-c("BB","BB-ED")
  Zind<-plyr::llply(.data = 1:M, .fun = function(x){
    if(TP==TRUE){
      PPSdata$Z[[x]][,indH0]>qnorm(1-alpha)}else{
        PPSdata$Z[[x]][,indH0]<qnorm(1-alpha)
      }
  }, .parallel = FALSE)
  
  FreqPPS<-plyr::llply(.data = 1:M, .fun = function(x){
    PPSdata$PPS[[x]][,indInt+indPPSH0][Zind[[x]]]
  }, .parallel = FALSE)
  
  ModelDF<-plyr::llply(.data = 1:M, .fun = function(x){
    rep(Models[x],sum(Zind[[x]]))
  }, .parallel = FALSE)
  
  DenPPS<-plyr::llply(.data = 1:M, .fun = function(x){
    density(FreqPPS[[x]], bw)
  }, .parallel = FALSE)
  
  DenPPSp<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(k){
    data.frame(x=DenPPS[[k]]$x[DenPPS[[k]]$x>0],
               y=DenPPS[[k]]$y[DenPPS[[k]]$x>0],
               variable=factor(rep(Models[k],sum(DenPPS[[k]]$x>0)),
                               levels=c("BB","BB-ED")))
  }, .parallel = FALSE))
  
  data2=data.frame(value=unlist(FreqPPS),
                   variable=factor(unlist(ModelDF),
                                   levels=c("BB","BB-ED")))
  
  phist <- gghistogram(
    data2, x = "value", 
    #add = "mean", rug = TRUE,
    fill = "variable",
    color= "variable", palette = colH,
    bins= 8,
    title = titlep,
    #xlab= expression(P(Z>z[alpha]~"|"~y[0]~","~y[1])),
    xlab=  "predicted probability of positive treatment effect",
    ylab = "frequency",
    alpha=.2
  )+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(legend.position = "none")
  
  # 2. Create the density plot with y-axis on the right
  # Remove x axis elements
  # pdensity <- ggdensity(
  #   data2, x = "value", 
  #   color= "variable", 
  #   palette = col,
  #   alpha = 0,
  #   lwd=1.5,
  #   linetype="variable"
  # ) +
  #   theme_half_open(11, rel_small = 1) +
  #   scale_linetype_manual(name="",values=linet, 
  #                         guide = guide_legend(override.aes=aes(fill=NA)))+
  #   #scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right")  +
  #   rremove("x.axis")+
  #   rremove("xlab") +
  #   rremove("x.text") +
  #   rremove("x.ticks") +
  #   rremove("y.axis")+
  #   rremove("ylab") +
  #   rremove("y.text") +
  #   rremove("y.ticks") +
  #   rremove("legend")
  
  pdensity <-ggplot(as.data.frame(DenPPSp),
                    aes(x=x,y=y,col=variable,linetype=variable))+geom_line(size=1.25)+
    scale_color_manual(values=col)+
    scale_linetype_manual(name="",values=linet, 
                          guide = guide_legend(override.aes=aes(fill=NA)))+
    theme_half_open(11, rel_small = 1) +
    rremove("x.axis")+
    rremove("xlab") +
    rremove("x.text") +
    rremove("x.ticks") +
    rremove("y.axis")+
    rremove("ylab") +
    rremove("y.text") +
    rremove("y.ticks") +
    rremove("legend")
  
  # 3. Align the two plots and then overlay them.
  aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
  ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])+
    theme(legend.position="none")
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

ROCfun<-function(pH0,pH1,cut=seq(0,1,length=1000)){
  x=apply(t(cut),2,function(x)mean(pH0>x,na.rm = T))
  y=apply(t(cut),2,function(x)mean(pH1>x,na.rm = T))
  out=data.frame(fpr=x,
                 tnr=y)
  return(out)
}

plot.ROC<-function(PPSdata,indInt=1,smoothind=T,alpha=.05,
                   col=c(customgreen,customgreen),
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
  
  print(round(dataAUC,3))
  
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
    theme(legend.position = "none")+
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


plot.false.stop <- function(PPSdata,bound=.1,futility=T,alpha=.05,
                            n.int,
                            col=c(customorange,customorange),
                            linet=c("solid","dashed"),
                            family="binomial",
                            method = "glm",
                            formula = 'y ~ x',
                            se=F){
  
  M<-length(PPSdata$PPS)
  
  data<-do.call("rbind",plyr::llply(.data = 1:M, .fun = function(x){
    if(futility==T){
      sapply(1:5,function(k) 
        mean(PPSdata$PPS[[x]][,k][PPSdata$Z[[x]][,1]>qnorm(1-alpha)]<bound))
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

formatTrial <- function(TrialsH0){
  if (colnames(TrialsH0)[1] == "X" ) {
    d = 1
  }else d = 0
  data <- vector("list", dim(TrialsH0)[2]/2)
  for (i in 0:(dim(TrialsH0)[2]/2-1)) {
    A <- as.data.frame(TrialsH0[,d+2*i+1:2])
    colnames(A) <- c("PPS.ID","PPS.ED")
    data[[i+1]] <- A
  }
  return(data)
}

