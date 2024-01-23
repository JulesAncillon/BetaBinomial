library(ggplot2)
library(readr)

sims = 2 # The number of in-sillico trials
Same_Pop = T # Can be True (T) or False (F)
interaction = F # With interactions (T) or without (F)
# Parameters trial
Pop_Type = "Diff" # Can be the Same ("Same") or Different ("Diff)
if (Same_Pop == T) {
  Pop_Type = "Same"
  P.X.ID = c(.5,.5,.5)
  P.X.ED = c(.5,.5,.5)
} else{
  P.X.ID = c(.7,.3,.6)
  P.X.ED = c(.2,.7,.7)
}
Interaction_Type = ""
if (interaction == T) {
  Interaction_Type = "_Int"
}


TrialsH0.C <- as.data.frame(matrix(0, ncol = 10, nrow = sims))
TrialsH1.C <- as.data.frame(matrix(1, ncol = 10, nrow = sims))
D.EC = read_csv(paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/Simulated/Data/",Pop_Type,Interaction_Type,"Ext_Data.csv",sep=""))

t1 <- Sys.time()
for (k in 1:1) {
  print(k)
  for (x in 1:sims) {
    D.Trials <- read_csv(paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/Simulated/Data/",Pop_Type,Interaction_Type,"Trial_",x,".csv",sep=""))
    D.aux = D.Trials[D.Trials$X.A!=2 & D.Trials$X.A!="E",]
    PPS <- pred.z.log(D.EC,D.aux,n.int[k], P.X.ID, S=10) # pred.z.[bb/log/nn]
    TrialsH0.C[x,2*k-1] <- PPS[1]
    TrialsH0.C[x,2*k] <- PPS[2]
    
    D.aux = D.Trials[D.Trials$X.A!=1 & D.Trials$X.A!="E",]
    D.aux$X.A[D.aux$X.A==2] = 1
    PPS_ <- pred.z.log(D.EC,D.aux,n.int[k], P.X.ID) # pred.z.[bb/log/nn]
    TrialsH1.C[x,2*k-1] <- PPS_[1]
    TrialsH1.C[x,2*k] <- PPS_[2]
  }
}
print(Sys.time()-t1)



####################
#Check if it looks nice
####################

ZH0 = read_csv(paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/Simulated/Data/",Pop_Type,Interaction_Type,"_ZH0.csv",sep=""))
ZH1 = read_csv(paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/Simulated/Data/",Pop_Type,Interaction_Type,"_ZH1.csv",sep=""))


SAME0 <- formatTrial(TrialsH0.C)
SAME1 <- formatTrial(TrialsH1.C)
SamePPS <- data.PPS(ZH0,ZH1,SAME0,SAME1,n.int)
BS_SAME <- data.BS(SamePPS,alpha=.05,n.int,interims=5)
plot.BS(BS_SAME,ylim=c(0,.2), col=c("black","black"), linet=c("solid","dashed"))

####################
#Gather the Results
####################

data0 <- as.data.frame(TrialsH0.C)
name0 = paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/Simulated/Results/",Pop_Type,Interaction_Type,"_TrialsH0_BB.csv",sep="")
write_csv(data0,path=name0)

data1<- as.data.frame(TrialsH1.C)
name1 = paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/Simulated/Results/",Pop_Type,Interaction_Type,"_TrialsH1_BB.csv",sep="")
write_csv(data1,path=name1)
