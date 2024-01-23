# Library 
rm(list=ls())
library(ggplot2)
library(readr)

Same_Pop = F
interaction = T
sims   = 200 # Number of Trials

# Parameters trial
Pop_Type = "Diff"
P.X.ED = c(.7,.3,.6) # Proba for the covariate of X.ED  c(.5,.5,.5) vs c(.6,.3,.5)  When interaction X3 ED = .8 and X3 ID = 5
P.X.ID = c(.2,.7,.7) # Proba for the covariate of X.ID vs c(.3,.6,.5)
if (Same_Pop == T) {
  Pop_Type = "Same"
  P.X.ED = c(.5,.5,.5) # Proba for the covariate of X.ED  c(.5,.5,.5) vs c(.6,.3,.5)
  P.X.ID = c(.5,.5,.5) # Proba for the covariate of X.ID vs c(.3,.6,.5)
  
}
b.IC   = c(qnorm(.7),-1,1,0)  # Beta*X.ID ; b.IC[1] is the intercept. 
b.EC   = c(qnorm(.7),-1,1,0) # Beta*X.ED ; b.EC[1] is the intercept. 
b.TE   = c(0,0,qnorm(.7)-qnorm(.5)) #Including Control arm  # Here 2nd Exp group has positive treatment effect, and 1st doesn't  ## Initially c(0,0,qnorm(.7)-qnorm(.5))
n.EC   = 2000 # Number of External observations 
n.C    = 350*1 # Number of observations for the control treamtment
n.E    = 350*2 # Number of observations for each experimental groups. Here random ration1:2:2
n.A    = length(b.TE) # Number of Arms
n.x    = length(P.X.ED) # Number of Biomarkers features
n.ID   = n.C+n.E*(n.A-1)  # Total number of 410 observations. 
n.int  = seq(50,250, by=100)
Interaction_Type = ""
if (interaction == T) {
  Interaction_Type = "_Int"
}
# Hyperparameters
a = 1
b = 1
alpha = .05

# Generate external and internal data
X.ED   = X.sample(n.EC,P.X.ED)  # Generate the covariates for the ED
P.Y.ED = pnorm(b.EC[1] + X.ED%*%b.EC[-1])  # Generate the probability for 1 or 0 in Y. 
Y.ED   = as.numeric(output(X.ED,c(0), b.EC, interaction))  #1*(runif(n.EC) <= P.Y.ED)  # Generate Y, according to the probability defined above
D.EC   = data.frame(arr.time=0,X.A="E",X=X.ED,ID=0,Y=Y.ED)
colnames(D.EC)<-c("arr.time","X.A","X1","X2","X3","ID","Y")
D.IC   = replicate(n=sims,sim.simple.trial.arm(c(n.C,n.E,n.E),P.X.ID,lambda=10,b.TE,b.IC, interaction), simplify=F) # 2 Experimental Arms


# Z-value
ZH0  <- do.call("rbind",plyr::llply(.data = 1:sims, .fun = function(x){
  z.test(D.IC[[x]]$Y[D.IC[[x]]$X.A!=2],D.IC[[x]]$X.A[D.IC[[x]]$X.A!=2])
}, .parallel = FALSE))
Z0 = mean(ZH0>qnorm(.95), na.rm = TRUE) 
ZH1 <- do.call("rbind",plyr::llply(.data = 1:sims, .fun = function(x){
  z.test(D.IC[[x]]$Y[D.IC[[x]]$X.A!=1],D.IC[[x]]$X.A[D.IC[[x]]$X.A!=1])
}, .parallel = FALSE))
Z1 = mean(ZH1>qnorm(.95), na.rm = TRUE)

print(paste("Z0 = ", Z0, " and Z1 = ",Z1))


# Register Data
for (i in 1:200) {
  trial <- paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/Simulated/Data/",Pop_Type,Interaction_Type,"_Trial_",i,".csv",sep="")
  write_csv(D.IC[[i]],path=trial)
}

Ext_Data <- paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/Simulated/Data/",Pop_Type,Interaction_Type,"Ext_Data.csv",sep="")
write_csv(D.EC,path=Ext_Data)
ZH0_ <- paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/Simulated/Data/",Pop_Type,Interaction_Type,"_ZH0.csv",sep="")
write_csv(as.data.frame(ZH0),path=ZH0_)
ZH1_ <- paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/Simulated/Data/",Pop_Type,Interaction_Type,"_ZH1.csv",sep="")
write_csv(as.data.frame(ZH1),path=ZH1_)


