# Neural Network
library(ANN2)
library(PoissonBinomial)

# Outcome.EC = D.EC
# x = 1
# D.aux = D.Trials[D.Trials$X.A!=2 & D.Trials$X.A!="E",]
# D.aux$X.A[D.aux$X.A==2] = 1
# Outcome.IC = D.aux
# NInt = n.int[1]

pred.z.nn <- function(Outcome.EC,Outcome.IC, NInt, alpha = 0.05, cont=F, L1 = 0, L2 = 0, P.X.ID = c(.5,.5,.5)){
  #Number of arms and total number of patients
  k = max(Outcome.IC$X.A)
  hiddenLayers = rep(5,5)                                # Neural network number of nodes 
  ## Interim ##
  ## Subsetting data to look ##
  Data.Interim =Outcome.IC[1:NInt,]
  n.Data      = data.frame(Arm=0:k,
                           #Interim
                           N=data.frame(table(Data.Interim$X.A))$Freq,             # A vector with the number of patients in each arms 
                           y=c(unlist(plyr::llply(.data = 0:k, .fun = function(x){ # A vector with the number of 1-outputs for each arms  #nb of success for each arms 
                             sum(Data.Interim$Y[Data.Interim$X.A==x])
                           }, .parallel = FALSE))),
                           #Final
                           N.fin=data.frame(table(Outcome.IC$X.A))$Freq,           # A vector with the number of patients in each arms 
                           y.fin=c(unlist(plyr::llply(.data = 0:k, .fun = function(x){ # The number of 1-outputs for each arms 
                             sum(Outcome.IC$Y[Outcome.IC$X.A==x])
                           }, .parallel = FALSE))))
  
  # Create the Training sets
  Train_C <- subset(Data.Interim[Data.Interim$X.A == 0,], select = -c(ID,X.A,Y.A0,Y.A1,Y.A2,arr.time))
  Train_C_ED <- as.data.frame(rbind(Train_C,subset(Outcome.EC, select = -c(ID,X.A,arr.time))))
  Train_E <- as.data.frame(subset(Data.Interim[Data.Interim$X.A != 0,], select = -c(ID,X.A,Y.A0,Y.A1,Y.A2,arr.time)))
  
  PPS_ID <- 0
  PPS_ED <- 0
  for (variable in vector) {
    # Prediction set 
    Pred_C <- as.data.frame(X.sample(n.Data$N.fin[1]-n.Data$N[1], P.X.ID,cont)) 
    Pred_E <- as.data.frame(X.sample(n.Data$N.fin[2]-n.Data$N[2], P.X.ID,cont)) 
    
    # Internal Data only Prediction
    NN_C = neuralnetwork(Train_C[,-4], Train_C$Y, hiddenLayers, loss.type = "log", L1 = 0,L2 = 0,batch.size = 10,
                         regression = FALSE, activ.functions = "sigmoid",n.epochs = 2000,learn.rates = 0.01)
    p_ID_C <- predict(NN_C, as.matrix(matrix(as.numeric(as.matrix(Pred_C)),dim(Pred_C)[1])))$probabilities[,2]
    
    NN_C_ED = neuralnetwork(Train_C_ED[,-4], Train_C_ED$Y, hiddenLayers, loss.type = "log", L1 = 0,L2 = 0,batch.size = 10,
                            regression = FALSE, activ.functions = "sigmoid",n.epochs = 2000,learn.rates = 0.01)
    p_ED_C <- predict(NN_C_ED, as.matrix(matrix(as.numeric(as.matrix(Pred_C)),dim(Pred_C)[1])))$probabilities[,2]
    
    NN_E = neuralnetwork(Train_E[,-4], Train_E$Y, hiddenLayers, loss.type = "log", L1 = 0,L2 = 0,batch.size = 10,
                         regression = FALSE, activ.functions = "sigmoid",n.epochs = 2000,learn.rates = 0.01)
    p_E <- predict(NN_E, as.matrix(matrix(as.numeric(as.matrix(Pred_E)),dim(Pred_E)[1])))$probabilities[,2]
    
    # Gives the values of the Test for every possible events in the 2 arms.
    zTestPred <- sapply(c(n.Data$y[2]:(n.Data$N.fin[2]-n.Data$N[2]+n.Data$y[2])), 
                        function(x) mapply(z.test2,c(n.Data$y[1]:(n.Data$N.fin[1]-n.Data$N[1]+n.Data$y[1])),x,n.Data$N.fin[1],n.Data$N.fin[2]))
    
    # Assuming success follows a Beta-binomial dist; this fct is looking at the joint probability of every possible events in the 2 arms. 
    PredProb.ID <- outer(dpbinom(x=NULL, as.numeric(unlist(p_ID_C)), method = "DivideFFT"), dpbinom(x=NULL, as.numeric(unlist(p_E)), method = "DivideFFT"), FUN = "*")
    PredProb.ED <- outer(dpbinom(x=NULL, as.numeric(unlist(p_ED_C)), method = "DivideFFT"), dpbinom(x=NULL, as.numeric(unlist(p_E)), method = "DivideFFT"), FUN = "*")
    
    id <- zTestPred>qnorm(1-alpha)
    PPS_ID <- PPS_ID + sum(PredProbC*id)
    PPS_ED <- PPS_ID + sum(PredProbE*id)
    print(PPS_ID)
  }
  PPS <- c(PPS_ID/S, PPS_ED/S)
  return(PPS) 
}
