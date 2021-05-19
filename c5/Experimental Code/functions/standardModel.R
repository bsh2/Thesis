standardLoss <- function(pars){
  squaredResiduals <- numeric(length(standardCriterion))
  inputData <- data.frame(days = 1:length(loads), 
                          loads = loads)
  for (i in 1:length(standardCriterion)){
    inputSubset <- inputData[1:standardIndex[i], ]
    modelledPerformance <- pars[1] + 
      pars[2]*(sum(inputSubset$loads * 
                   exp(-(standardIndex[i] - inputSubset$days) / pars[3]))) -
      pars[4]*(sum(inputSubset$loads * 
                   exp(-(standardIndex[i] - inputSubset$days) / pars[5])))
    squaredResiduals[i] <- (modelledPerformance - standardCriterion[i])^2  
  }
  return(sum(squaredResiduals))
}

standardFit <- function(vals){
  startVals <- as.numeric(vals[1:5])
  startRSS <- as.numeric(vals[6])
  parsFound <- optim(par = startVals,
                     fn = standardLoss,
                     method = "L-BFGS-B",
                     lower = standardConstraints$lower,
                     upper = standardConstraints$upper,
                     hessian = TRUE,
                     control = list(
                       maxit = 10000,
                       trace = FALSE
                     ))
  tempPars <- as.numeric(parsFound$par)
  tempParDiff <- tempPars - standardPars
  tempLoss <- as.numeric(parsFound$value)
  tempCounts <- as.numeric(parsFound$counts)
  tempConvergence <- as.numeric(parsFound$convergence)
  tempEigen <- eigen(parsFound$hessian)
  tempEigen <- as.numeric(tempEigen$values)
  tempPredictions <- standardCompute(tempPars, loads, 
                                     returnObject = "performance")
  tempPredictions <- tempPredictions[standardIndex]
  tempRsq <- RSQfunc(tempPredictions, standardCriterion)
  tempRMSE <- RMSEfunc(tempPredictions, standardCriterion, 
                       length(standardIndex))
  tempMAPE <- MAPEfunc(tempPredictions, standardCriterion)
  returnVector <- c(startVals, startRSS, tempPars, tempLoss, tempCounts, 
                    tempConvergence, tempEigen, tempRsq, 
                    tempRMSE, tempMAPE, tempParDiff)
  return(returnVector)
}

standardCompute <- function(pars, loads, returnObject = "all"){
  pars <- as.numeric(pars)
  p <- numeric(length = length(loads))
  fitness <- numeric(length = length(loads))
  fatigue <- numeric(length = length(loads))
  s <- 1:length(loads)
  df0 <- data.frame(s, "ws" = loads)
  for (n in 1:length(s)){
    df1 <- df0[1:s[n], ]
    fitness[n] <- pars[2] * sum( df1$ws * exp(- (n - df1$s) / pars[3]) )
    fatigue[n] <- pars[4] * sum( df1$ws * exp(- (n - df1$s) / pars[5]) )
    p[n] <- pars[1] + fitness[n] - fatigue[n]
  }
  if (returnObject == "performance"){
    return(p)} else{
      return(data.frame("fitness" = fitness, "fatigue" = fatigue,
                        "performance" = p))  
    }
}