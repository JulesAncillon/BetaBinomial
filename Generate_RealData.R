# Library 
library(extraDistr)
library(utils)
library(dplyr)
library(grid)
library(ggplot2)
library(caret)
library(varhandle)
library(aod)
library(ranger)
library(tidyverse)

# Source
source("/Users/julesancillon/Desktop/Harvard/MastersThesis/Pred/Predict.fun.R")

# Read Data
d1 <- read_csv("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/GBM/GBM_Data/AugmentedCORE.csv")
d2 <- read_csv("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/GBM/GBM_Data/AugmentedCENTRIC.csv")

core <- select(.data = d1, c("age","sex","kps","mgmt","eor","os_status"))
centric <- select(.data = d2, c("age","sex","kps","mgmt","eor","os_status"))
core <- core[complete.cases(core), ]
centric <- centric[!(centric$kps=="50-60"),]
centric <- centric[complete.cases(centric), ]

# Select Parameters
Same_Pop = F
sims = 200
n.int  = seq(100,500, by=100)

# Data Generation
Pop_Mix = .1
Pop_Type = "Diff"
P.TE = 0.70
if (Same_Pop == T) {
  Pop_Mix = .5
  Pop_Type = "Same"
  P.TE = 0.55
}
A1 <- core[1:(550*(1-Pop_Mix)),]
A2 <- core[(550*(1-Pop_Mix)+1):550,]
B1 <- centric[1:(550*(1-Pop_Mix)),]
B2 <- centric[(550*(1-Pop_Mix)+1):550,]
C1 <- rbind(A1,B2) %>% drop_na()
C2 <- rbind(A2,B1) %>% drop_na()

Datas <- list(C1,C2)
N = dim(C1)[1]
N0 = trunc(dim(C1)[1]*.66,0)

D.Trials = replicate(n=sims,simulate.trial(Datas, N, N0,P.TE), simplify=F) 

# Z-value
ZH0  <- do.call("rbind",plyr::llply(.data = 1:sims, .fun = function(x){
  z.test(D.Trials[[x]]$Y[D.Trials[[x]]$X.A!=2 & D.Trials[[x]]$X.A!="E"],D.Trials[[x]]$X.A[D.Trials[[x]]$X.A!=2 & D.Trials[[x]]$X.A!="E"])
}, .parallel = FALSE))
Z0 = mean(ZH0>qnorm(.95), na.rm = TRUE) 
ZH1  <- do.call("rbind",plyr::llply(.data = 1:sims, .fun = function(x){
  z.test(D.Trials[[x]]$Y[D.Trials[[x]]$X.A!=1 & D.Trials[[x]]$X.A!="E"],D.Trials[[x]]$X.A[D.Trials[[x]]$X.A!=1 & D.Trials[[x]]$X.A!="E"])
}, .parallel = FALSE))
Z1 = mean(ZH1>qnorm(.95), na.rm = TRUE)
print(paste("Z0 = ", Z0, " and Z1 = ",Z1))


# Register Data
for (i in 1:sims) {
  trial <- paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/Real/Data/",Pop_Type,"Trial_",i,sep="")
  write_csv(D.Trials[[i]],path=trial)
}

ZH0_ <- paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/Real/Data/",Pop_Type,"_ZH0.csv")
write_csv(as.data.frame(ZH0),path=ZH0_)
ZH1_ <- paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/Real/Data/",Pop_Type,"_ZH1.csv")
write_csv(as.data.frame(ZH1),path=ZH1_)


