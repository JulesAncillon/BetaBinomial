library(PoissonBinomial)
# Outcome.EC = D.EC
# x = 1
# Outcome.IC = D.IC[[x]][D.IC[[x]]$X.A!=2,]
# NInt = n.int[3]

pred.z.log <- function(Outcome.EC,Outcome.IC, NInt, alpha = 0.05,cont = F, P.X.ID = c(.4,.7,.5), S = 100){
  #Number of arms and total number of patients
  k = max(Outcome.IC$X.A)
  
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
  
  # Training set TrialsH0
  Train_C <- as.data.frame(subset(Data.Interim[Data.Interim$X.A == 0,], select = -c(X.A,ID,Y.A0,Y.A1,Y.A2,arr.time)))
  Train_C_ED <- as.data.frame(rbind(Train_C,subset(Outcome.EC, select = -c(X.A,ID,arr.time))))
  Train_E <- as.data.frame(subset(Data.Interim[Data.Interim$X.A != 0,], select = -c(X.A,ID,Y.A0,Y.A1,Y.A2,arr.time)))
  
  PPS_ID <- 0
  PPS_ED <- 0
  for (s in 1:S) {
    # Prediction set 
    Pred_C <- as.data.frame(X.sample(n.Data$N.fin[1]-n.Data$N[1], P.X.ID, cont)) 
    Pred_E <- as.data.frame(X.sample(n.Data$N.fin[2]-n.Data$N[2], P.X.ID, cont)) 
    
    # Creating the model TrialsH0
    Log_C <- glm(Train_C$Y ~ ., data = Train_C, family = "binomial")
    Log_C_ED <- glm(Train_C_ED$Y ~ ., data = Train_C_ED, family = "binomial")
    Log_E <- glm(Train_E$Y ~ ., data = Train_E, family = "binomial")
    
    # Using the model to find the probabilities 
    pC <- predict(Log_C, newdata = Pred_C, type = "response")
    pC_ED <- predict(Log_C_ED, newdata = Pred_C, type = "response")
    pE<- predict(Log_E, newdata = Pred_E, type = "response")
    
    # Gives the values of the Test for every possible events in the 2 arms.
    zTestPred <- sapply(c(n.Data$y[2]:(n.Data$N.fin[2]-n.Data$N[2]+n.Data$y[2])), 
                        function(x) mapply(z.test2,c(n.Data$y[1]:(n.Data$N.fin[1]-n.Data$N[1]+n.Data$y[1])),x,n.Data$N.fin[1],n.Data$N.fin[2]))
    
    
    # Assuming success follows a Beta-binomial dist; this fct is looking at the joint probability of every possible events in the 2 arms. 
    PredProbC <- outer(dpbinom(x=NULL, pC, method = "DivideFFT"), dpbinom(x=NULL, pE, method = "DivideFFT"), FUN = "*")
    PredProbE <- outer(dpbinom(x=NULL, pC_ED, method = "DivideFFT"), dpbinom(x=NULL, pE, method = "DivideFFT"), FUN = "*")
    
    id <- zTestPred>qnorm(1-alpha)
    PPS_ID <- PPS_ID + sum(PredProbC*id)
    PPS_ED <- PPS_ID + sum(PredProbE*id)
    print(PPS_ID)
  }
  
  PPS <- c(PPS_ID/S, PPS_ED/S)
  return(PPS)
}



