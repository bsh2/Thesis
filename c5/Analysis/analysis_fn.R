standardLoss <- function(pars, standardCriterion, loads, standardIndex){
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
    return(p)} 
  if (returnObject == "fitness"){
    return(fitness)
  }
  if (returnObject == "fatigue"){
    return(fatigue)
  }
  if (returnObject == "all"){
    return(data.frame("fitness" = fitness, "fatigue" = fatigue,
                      "performance" = p))  
  }
}

standard_true_pars <- c(100, 0.72, 28.5, 1.2, 8.6)
loads <- read.csv("loads.csv")
loads <- loads$load
standardCriterion1 <- standardCompute(standard_true_pars, loads, returnObject = "performance")
standardIndex1 <- seq(1, length(standardCriterion1), by = 1)
standardCriterion2 <- standardCriterion1[seq(1,length(standardCriterion1),by=2)]
standardIndex2 <- seq(1, length(standardCriterion1), by = 2)
standardCriterion3 <- standardCriterion1[seq(1,length(standardCriterion1),by=3)]
standardIndex3 <- seq(1, length(standardCriterion1), by = 3)

calvertLoss <- function(pars, calvertCriterion, loads, calvertIndex){
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
    return(p)} 
  if (returnObject == "fitness"){
    return(fitness)
  }
  if (returnObject == "fatigue"){
    return(fatigue)
  }
  if (returnObject == "all"){
    return(data.frame("fitness" = fitness, "fatigue" = fatigue,
                      "performance" = p))  
  }
}

calvert_true_pars <- c(100, 0.72, 32.5, 4.3, 1.05, 8.6)
calvertCriterion1 <- calvertCompute(calvert_true_pars, loads, returnObject = "performance")
calvertIndex1 <- seq(1, length(calvertCriterion1), by = 1)
calvertCriterion2 <- calvertCriterion1[seq(1,length(calvertCriterion1),by=2)]
calvertIndex2 <- seq(1, length(calvertCriterion1), by = 2)
calvertCriterion3 <- calvertCriterion1[seq(1,length(calvertCriterion1),by=3)]
calvertIndex3 <- seq(1, length(calvertCriterion1), by = 3)
