#' Simulate platform hybrid trial that adds arms, stops for futility and incorporate external data
#' @param Outcome.IC Possible outcomes of internal trial (including arrival times, patient characteristics)
#' @param Outcome.EC Outcomes of Observed external trial
#' @param NInt Number of patients at interim
#' @param a Hyperparameter BB
#' @param b Hyperparameter BB
#' @export
#' @examples
#' pred.z.bb()
pred.z.bb = function(Outcome.EC,Outcome.IC,NInt,a=1,b=1,alpha=.05){
  
  #Number of arms and total number of patients
  k = max(Outcome.IC$X.A)
  ## Interim ##
  ## Subsetting data to look ##
  Data.Interim=Outcome.IC[1:NInt,]
  
  # summary - nr per arm (N, status)
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
                           }, .parallel = FALSE))))#,
  #Status=rep(0,(k+1)))
  
  # Gives the values of the Test for every possible events in the 2 arms.
  zTestPred <- sapply(c(n.Data$y[2]:(n.Data$N.fin[2]-n.Data$N[2]+n.Data$y[2])), 
                      function(x) mapply(z.test2,c(n.Data$y[1]:(n.Data$N.fin[1]-n.Data$N[1]+n.Data$y[1])),x,n.Data$N.fin[1],n.Data$N.fin[2]))
  
  # Assuming success follows a Beta-binomial dist; this fct is looking at the joint probability of every possible events in the 2 arms. 
  # PredProb <- outer(extraDistr::dbbinom(0:(n.Data$N.fin[1]-n.Data$N[1]),(n.Data$N.fin[1]-n.Data$N[1]),a+n.Data$y[1],b+n.Data$N[1]-n.Data$y[1]),
  #                   extraDistr::dbbinom(0:(n.Data$N.fin[2]-n.Data$N[2]),(n.Data$N.fin[2]-n.Data$N[2]),a+n.Data$y[2],b+n.Data$N[2]-n.Data$y[2]),
  #                   FUN = "*")
  
  #sum(PredProb) #1 Check
  PredProbE <- outer(extraDistr::dbbinom(0:(n.Data$N.fin[1]-n.Data$N[1]),(n.Data$N.fin[1]-n.Data$N[1]),a+n.Data$y[1]+sum(Outcome.EC$Y),b+n.Data$N[1]+nrow(Outcome.EC)-n.Data$y[1]-sum(Outcome.EC$Y)),
                     extraDistr::dbbinom(0:(n.Data$N.fin[2]-n.Data$N[2]),(n.Data$N.fin[2]-n.Data$N[2]),a+n.Data$y[2],b+n.Data$N[2]-n.Data$y[2]),
                     FUN = "*")
 
   #sum(PredProbE) #1
  id <- zTestPred>qnorm(1-alpha)
  
  PPS <- data.frame(PPS.ED = sum(PredProbE*id)) #, PPS.ID = sum(PredProb*id))
                    
  if(sum(is.na(PPS))!=0){
    # PPS$PPS.ID = .5 
    PPS$PPS.ED = .5
    }

  return(PPS) 
}
