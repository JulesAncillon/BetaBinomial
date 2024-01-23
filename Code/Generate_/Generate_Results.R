library(readr)
###############
# Parameters 
###############
Same_Pop = 'D'  # D for Different and S for Same 
Data_Type = 'R' # S for Simulated and R for Real 
Interactions = T

# Machine Gun 
Interaction_Type = ""
Pop_Type = "Diff"
if (Same_Pop == "S") {
  Pop_Type = "Same"
}
if (Interactions == T & Data_Type == "S") {
  Interaction_Type = "_Int"
}
Type = "Real"
n.int = seq(100,500, by=100)
if (Data_Type == 'S') {
  Type = "Simulated"
  n.int = seq(50,250, by=50)
}

print(paste(Type," Data with ",Pop_Type," biomarkers", "and Interactions : ",Interactions ))


###############
# Read Data
###############
ZH0 = read_csv(paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/",Type,"/Data/",Pop_Type,Interaction_Type,"_ZH0.csv",sep=""))
ZH1 = read_csv(paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/",Type,"/Data/",Pop_Type,Interaction_Type,"_ZH1.csv",sep=""))

TrialsH0_BB = formatTrial(read_csv(paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/",Type,"/Results/",Pop_Type,Interaction_Type,"_TrialsH0_BB.csv",sep="")))
TrialsH1_BB = formatTrial(read_csv(paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/",Type,"/Results/",Pop_Type,Interaction_Type,"_TrialsH1_BB.csv",sep="")))
TrialsH0_LOG = formatTrial(read_csv(paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/",Type,"/Results/",Pop_Type,Interaction_Type,"_TrialsH0_LOG.csv",sep="")))
TrialsH1_LOG = formatTrial(read_csv(paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/",Type,"/Results/",Pop_Type,Interaction_Type,"_TrialsH1_LOG.csv",sep="")))
TrialsH0_NN1 = formatTrial(read_csv(paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/",Type,"/Results/",Pop_Type,Interaction_Type,"_TrialsH0_NN1.csv",sep="")))
TrialsH1_NN1 = formatTrial(read_csv(paste("/Users/julesancillon/Desktop/Harvard/MastersThesis/Data/Final_Paper/",Type,"/Results/",Pop_Type,Interaction_Type,"_TrialsH1_NN1.csv",sep="")))

PPS_BB = data.PPS(ZH0,ZH1,TrialsH0_BB,TrialsH1_BB,n.int)
PPS_LOG = data.PPS(ZH0,ZH1,TrialsH0_LOG,TrialsH1_LOG,n.int)
PPS_NN = data.PPS(ZH0,ZH1,TrialsH0_NN1,TrialsH1_NN1,n.int)

###############
# Brier Scores
###############

# Single
PPS <- data.PPS(ZH0,ZH1,TrialsH0_LOG,TrialsH1_LOG,n.int)
BS <- data.BS(PPS,alpha=.05,n.int,interims=5)
#plot.BS(BS,ylim=c(0,.5), col=c('#999999','#999999'), linet=c("solid","dashed"))

# Triple
dataBS <- data.BS_triple(PPS_BB,PPS_LOG,PPS_NN,alpha=.05,n.int,interims=5)
plot.BS_triple(dataBS)


###############
# False Stop
###############
bound=.1

plot.false.stop_triple(PPS_BB,PPS_LOG,PPS_NN,bound=bound,futility=T,alpha=.05,n.int)
#
#plot.false.stop(PPS_NN,bound=.1,futility=T,alpha=.05,n.int)
#plot.stop(PPS_LOG,futility=T,alpha=.05,n.int)

###############
# True Stop
###############
bound=.1

plot.true.stop_triple(PPS_BB,PPS_LOG,PPS_NN,bound=bound,futility=T,alpha=.05,n.int)
