calvertLoss <- function(pars){
  squaredResiduals <- numeric(length(calvertCriterion))
  inputData <- data.frame(days = 1:length(loads), 
                          loads = loads)
  for (i in 1:length(calvertCriterion)){
    inputSubset <- inputData[1:calvertIndex[i], ]
    modelledPerformance <- pars[1] + 
      pars[2]*(sum(inputSubset$loads * 
                     (exp(-(calvertIndex[i] - inputSubset$days) / pars[3]) -
                      exp(-(calvertIndex[i] - inputSubset$days) / pars[4])))) -
      pars[5]*(sum(inputSubset$loads * 
                     exp(-(calvertIndex[i] - inputSubset$days) / pars[6])))
    squaredResiduals[i] <- (modelledPerformance - calvertCriterion[i])^2  
  }
  return(sum(squaredResiduals))
}

calvertFit <- function(vals){
  startVals <- as.numeric(vals[1:6])
  startRSS <- as.numeric(vals[7])
  parsFound <- optim(par = startVals,
                     fn = calvertLoss,
                     method = "L-BFGS-B",
                     lower = calvertConstraints$lower,
                     upper = calvertConstraints$upper,
                     hessian = TRUE,
                     control = list(
                       maxit = 10000,
                       trace = FALSE
                     ))
  tempPars <- as.numeric(parsFound$par)
  tempParDiff <- tempPars - calvertPars
  tempLoss <- as.numeric(parsFound$value)
  tempCounts <- as.numeric(parsFound$counts)
  tempConvergence <- as.numeric(parsFound$convergence)
  tempEigen <- eigen(parsFound$hessian)
  tempEigen <- as.numeric(tempEigen$values)
  tempPredictions <- calvertCompute(tempPars, loads, 
                                    returnObject = "performance")
  tempPredictions <- tempPredictions[calvertIndex]
  tempRsq <- RSQfunc(tempPredictions, calvertCriterion)
  tempRMSE <- RMSEfunc(tempPredictions, calvertCriterion, 
                       length(calvertIndex))
  tempMAPE <- MAPEfunc(tempPredictions, calvertCriterion)
  returnVector <- c(startVals, startRSS, tempPars, tempLoss, tempCounts, 
                    tempConvergence, tempEigen, tempRsq, 
                    tempRMSE, tempMAPE, tempParDiff)
  return(returnVector)
}

calvertCompute <- function(pars, loads, returnObject = "all"){
  pars <- as.numeric(pars)
  p <- numeric(length = length(loads))
  fitness <- numeric(length = length(loads))
  fatigue <- numeric(length = length(loads))
  s <- 1:length(loads)
  df0 <- data.frame(s, "ws" = loads)
  for (n in 1:length(s)){
    df1 <- df0[1:s[n], ]
    fitness[n] <- pars[2] * sum( df1$ws * (exp(- (n - df1$s) / pars[3])- 
                                           exp(- (n - df1$s) / pars[4])))
    fatigue[n] <- pars[5] * sum( df1$ws * exp(- (n - df1$s) / pars[6]) )
    p[n] <- pars[1] + fitness[n] - fatigue[n]
  }
  if (returnObject == "performance"){
    return(p)} else{
      return(data.frame("fitness" = fitness, "fatigue" = fatigue,
                        "performance" = p))  
    }
}