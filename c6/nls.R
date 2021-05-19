# ******************************************************************************
# Language: R (r-project.org.uk)
# Author: Ben Stephens Hemingway
# License: GNU GPL v3
# Thesis chapter: 6
# Subsection: 6.2.1
# Description: Nonlinear least-squares (standard model, fitness-delay, VDR)
# R version 4.04 (2021-02-15 "Lost Library Book")
# Package dependencies:
library(optimx)            # Standard optimization package
# ******************************************************************************

# ------------------------------------------------------------------------------
# 6.1 | Objective function (RSS): Standard model - For loop approach
# ------------------------------------------------------------------------------

standardObjectiveSS <- function(pars, loads, perfVals, initial = FALSE, 
                                maximise = FALSE){
  
  # (pars) NOTES:
  # ----------------------------------------------------------------------------
  # Pos.     [1] [2] [3] [4]  [5] [6] [7]
  # pars = c(p*, kg, Tg, kh, Th)          (if initial = FALSE)
  # pars = c(p*, kg, Tg, kh, Th, qg, qh)  (if initial = TRUE)
  # ----------------------------------------------------------------------------
  
  nMeasurements <- length(perfVals$performance) # Number of performance measurements
  
  # Zeroed vector of length equal to number of performance measurements
  squaredResiduals <- numeric(length = nMeasurements)
  
  # For each performance measurement calculate (modeled - measured)^2 under pars
  for (n in 1:nMeasurements){
    
    dayT <- perfVals$day[n]              # Day of measured performance
    measured <- perfVals$performance[n]  # Measured performance value on dayT
    
    # Isolate the required load data to compute the model up to dayT
    # Note: 1:dayT rather than 1:(dayT - 1) as the first row in the loads array is w(0)=0
    inputSubset <- loads[1:dayT, ]
    
    if (initial == TRUE){        
      # Include residual effects (initial components)
      initFitness <- pars[6] * exp(-(dayT) / pars[3])
      initFatigue <- pars[7] * exp(-(dayT) / pars[5])
    } else{
      initFitness <- 0
      initFatigue <- 0
    }
    
    # Compute modelled performance on dayT
    model <- pars[1] + initFitness - initFatigue +
      pars[2] * (sum(inputSubset$load * exp(-(dayT - inputSubset$day) / pars[3]))) -
      pars[4] * (sum(inputSubset$load * exp(-(dayT - inputSubset$day) / pars[5])))
    
    # Compute the squared residual value (model - measured)^2
    squaredResiduals[n] <- (model - measured)^2  
    
  } # Loop updates until iterated over all available measurements
  
  # Output
  if(maximise == FALSE){
    return(sum(squaredResiduals))
  }
  if(maximise == TRUE){ 
    # For optimisation algorithms that maximise by default
    return(-1 * sum(squaredResiduals))
  }
}

# ------------------------------------------------------------------------------
# 6.2 | Modification of listing 6.1 to incorporate the fitness-delay model
# ------------------------------------------------------------------------------

fitnessDelayObjectiveSS <- function(pars, loads, perfVals, initial = FALSE,
                                    maximise = FALSE){
  
  # INPUT NOTES:
  # ---------------------------------------------------------
  #         [1] [2] [3] [4]  [5] [6] [7] [8]
  # Pars: c(p*, kg, Tg, Tg2, kh, Th)          initial = FALSE
  # Pars: c(p*, kg, Tg, Tg2, kh, Th, qg, qh)  initial = TRUE
  # ----------------------------------------------------------
  
  nMeasurements <- length(perfVals$performance) # Number of performance measurements
  
  # Zeroed vector of length equal to number of performance measurements
  squaredResiduals <- numeric(length = nMeasurements)
  
  # For each performance measurement calculate (modeled - measured)^2 under pars
  for (n in 1:nMeasurements){
    
    dayT <- perfVals$day[n]              # Day of measured performance
    measured <- perfVals$performance[n]  # Measured performance value on dayT
    
    # Isolate the required load data to compute the model up to dayT
    # Note: 1:dayT rather than 1:(dayT - 1) as the first row in the loads array is w(0)=0
    inputSubset <- loads[1:dayT, ]
    
    if (initial == TRUE){
      initFitness <- pars[7] * exp(-(dayT) / pars[3])
      initFatigue <- pars[8] * exp(-(dayT) / pars[6])
    } else{
      initFitness <- 0
      initFatigue <- 0
    }
    
    # Compute modelled performance on dayT under pars // p\hat(dayT) = p* + g(dayT) - h(dayT)
    model <- pars[1] + initFitness - initFatigue +
      pars[2] * (sum(inputSubset$load * (exp(-(dayT - inputSubset$day) / pars[3]) -
                                           exp(-(dayT - inputSubset$day) / pars[4])))) - 
      pars[5] * (sum(inputSubset$load * exp(-(dayT - inputSubset$day) / pars[6])))
    
    
    # Compute the squared residual value (model - measured)^2
    squaredResiduals[n] <- (model - measured)^2  
    
  } # Loop updates until iterated over all available measurements
  
  # Output
  if(maximise == FALSE){
    return(sum(squaredResiduals))
  }
  if(maximise == TRUE){
    return(-1 * sum(squaredResiduals))
  }
}

# ------------------------------------------------------------------------------
# 6.3 | Simulation function: Standard Model - For loop approach
# ------------------------------------------------------------------------------

simulateStandard <- function(pars, loads, initialPars = c(0,0), returnObject = "all"){
  
  # Set up zeroed vectors of required length
  seriesLength <- tail(loads$day, 1)
  performance <- numeric(length = seriesLength)     # Model performance
  fitness <- numeric(length = seriesLength)         # Model fitness
  fatigue <- numeric(length = seriesLength)         # Model fatigue
  initialFitness <- numeric(length = seriesLength)  # Residual effects (initial component)
  initialFatigue <- numeric(length = seriesLength)
  
  # Calculate model fitness g(t), fatigue h(t), and 
  # performance  p(t) for t = 1:seriesLength
  for (t in 1:seriesLength){
    
    # Isolate the required load data for calculating p(t) 
    # (i.e., loads from day 0 to day t-1)
    inputSubset <- loads[loads$day < t, ]
    
    # Residual effects from initial components at time point t
    initialFitness[t] <- initialPars[1] * exp(-(t) / pars[3])
    initialFatigue[t] <- initialPars[2] * exp(-(t) / pars[4])
    
    # Compute g(t), h(t), p(t) for current t
    fitness[t] <- pars[2] * sum(inputSubset$load * exp(- (t - inputSubset$day)/ pars[3]))
    fatigue[t] <- pars[4] * sum(inputSubset$load * exp(- (t - inputSubset$day)/ pars[5]))
    performance[t] <- pars[1] + fitness[t] - fatigue[t] + initialFitness[t] - initialFatigue[t]
    
  } # Loop index updates (t <- t+1) until t = seriesLength
  
  # Output
  if (returnObject == "performance"){return(performance)}
  if (returnObject == "fitness"){return(fitness)}
  if (returnObject == "fatigue"){return(fatigue)}
  if (returnObject == "all"){
    return(data.frame("day" = 1:seriesLength,
                      "initial_fitness" = initialFitness,
                      "initial_fatigue" = initialFatigue,
                      "fitness" = fitness, "fatigue" = fatigue,
                      "performance" = performance))
  }
  
} # End function (closing bracket)


# ------------------------------------------------------------------------------
# 6.4 | Simulation function: Fitness-delay model - For loop approach
# ------------------------------------------------------------------------------

simulateFitnessDelay <- function(pars, loads, initialPars = c(0,0), returnObject = "all"){
  
  # Set up zeroed vectors of required length
  seriesLength <- tail(loads$day, 1)
  performance <- numeric(length = seriesLength)     # Model performance
  fitness <- numeric(length = seriesLength)         # Model fitness
  fatigue <- numeric(length = seriesLength)         # Model fatigue
  initialFitness <- numeric(length = seriesLength)  # Residual effects (initial component)
  initialFatigue <- numeric(length = seriesLength)
  
  # Calculate model fitness g(t), fatigue h(t), and 
  # performance  p(t) for t = 1:seriesLength
  for (t in 1:seriesLength){
    
    # Isolate the required load data for calculating p(t) 
    # (i.e., loads from day 0 to day t-1)
    inputSubset <- loads[loads$day < t, ]
    
    # Residual effects from initial components at time point t
    initialFitness[t] <- initialPars[1] * exp(-(t) / pars[3])
    initialFatigue[t] <- initialPars[2] * exp(-(t) / pars[6])
    
    # Compute g(t), h(t), p(t) for current t
    fitness[t] <- pars[2] * sum(inputSubset$load * (exp(-(t - inputSubset$day) / pars[3])  - 
                                                    exp(-(t - inputSubset$day) / pars[4])))
    fatigue[t] <- pars[5] * sum(inputSubset$load * exp(-(t - inputSubset$day) / pars[6]))
    performance[t] <- pars[1] + fitness[t] - fatigue[t] + initialFitness[t] - initialFatigue[t]
    
  } # Loop index updates (t <- t+1) until t = seriesLength
  
  # Output
  if (returnObject == "performance"){return(performance)}
  if (returnObject == "fitness"){return(fitness)}
  if (returnObject == "fatigue"){return(fatigue)}
  if (returnObject == "all"){
    return(data.frame("day" = 1:seriesLength,
                      "initial_fitness" = initialFitness,
                      "initial_fatigue" = initialFatigue,
                      "fitness" = fitness, "fatigue" = fatigue,
                      "performance" = performance))
  }
  
} # End function (closing bracket)



# ------------------------------------------------------------------------------
# 6.5 | Objective function (RSS): VDR model - For loop approach
# ------------------------------------------------------------------------------

vdrObjectiveSS <- function(pars, loads, perfVals, initial = FALSE, maximise = FALSE){
  
  # INPUT NOTES:
  # ----------------------------------------------------------
  #         [1] [2] [3] [4] [5] [6]  [7] [8]
  # Pars: c(p*, kg, Tg, kh, Th, Th2)          initial = FALSE
  # Pars: c(p*, kg, Tg, kh, Th, Th2, qg, qh)  initial = TRUE
  # ----------------------------------------------------------
  
  nMeasurements <- length(perfVals$performance)  # Number of performance measurements
  
  # Zeroed vector of length equal to number of performance measurements
  squaredResiduals <- numeric(length = nMeasurements)
  
  # For each performance measurement calculate (modelled - measured)^2 under pars
  for (n in 1:nMeasurements){
    
    dayT <- perfVals$day[n]                # Day of measured performance
    measured <- perfVals$performance[n]    # Measured performance value on dayT
    
    # Isolate the required load data to compute the model up to dayT (i.e., from t=0 to dayT - 1)
    inputSubset <- loads[1:dayT, ]
    
    # Initial components
    if (initial == TRUE){
      initFitness <- pars[7] * exp(-(dayT) / pars[3])
      initFatigue <- pars[8] * exp(-(dayT) / pars[6])
    } else{
      initFitness <- 0
      initFatigue <- 0
    }
    
    # Set up a zeroed vector to hold the variable gain term values kh2(i) for i=0 to dayT - 1
    kh2 <- numeric(length = dayT)  # Variable gain term vector
    
    # Calculate the variable gain term kh2(i) for i=0,1,2,...,dayT-1 (Recursive)
    for (i in 1:dayT){
      kh2[i] <- sum( inputSubset$load[1:i] * exp( -((inputSubset$day[i]-inputSubset$day[1:i])/pars[6])))
    }
    
    # Compute modelled performance on dayT under pars // p\hat(dayT) = p* + g(dayT) - h(dayT)
    model <- pars[1] + initFitness - initFatigue +
      pars[2] * (sum(inputSubset$load * exp(-(dayT - inputSubset$day) / pars[3]))) - 
      pars[4] * (sum(kh2 * inputSubset$load * exp(-(dayT - inputSubset$day) / pars[5])))
    
    # Compute the squared residual value (model - measured)^2
    squaredResiduals[n] <- (model - measured)^2
    
  }
  
  # Output
  if(maximise == FALSE){
    return(sum(squaredResiduals))
  }
  if(maximise == TRUE){
    return(-1 * sum(squaredResiduals))
  }
}

# ------------------------------------------------------------------------------
# 6.7 | Simulation function: Fitness-delay model - For loop approach
# ------------------------------------------------------------------------------

simulateVDR <- function(pars, loads, initialPars = c(0,0), returnObject = "all"){
  
  # pars = c(p*, k_g, Tau_g, k_h, Tau_h, Tau_h2)
  
  # Set up zeroed vectors of required length
  seriesLength <- tail(loads$day, 1)
  performance <- numeric(length = seriesLength)     # Model performance
  fitness <- numeric(length = seriesLength)         # Model fitness
  fatigue <- numeric(length = seriesLength)         # Model fatigue
  initialFitness <- numeric(length = seriesLength)  # Residual effects (initial component)
  initialFatigue <- numeric(length = seriesLength)
  kh2Dat <- matrix(data = NaN, nrow = tail(loads$day, 1), ncol = tail(loads$day, 1))
  
  # Calculate model fitness g(t), fatigue h(t), and 
  # performance  p(t) for t = 1:seriesLength
  for (t in 1:seriesLength){
    
    # Isolate the required load data for calculating p(t) 
    # (i.e., loads from day 0 to day t-1)
    inputSubset <- loads[loads$day < t, ]
    
    # Residual effects from initial components at time point t
    initialFitness[t] <- initialPars[1] * exp(-(t) / pars[3])
    initialFatigue[t] <- initialPars[2] * exp(-(t) / pars[4])
    
    # Set up a zeroed vector to hold the variable gain term values kh2(i) for i=0 to dayT - 1
    kh2 <- numeric(length = t)  # Variable gain term vector
    
    # Calculate the variable gain term kh2(i) for i=0,1,2,...,dayT-1 (Recursive)
    for (i in 1:t){
      kh2[i] <- sum(inputSubset$load[1:i] * 
                    exp( -((inputSubset$day[i]-inputSubset$day[1:i])/pars[6])))
    }
    
    # For each iteration t+1 of the outer loop, save kh2(i) values where i = 0 to t-1
    kh2Dat[1:t, t] <- kh2
    
    # Compute g(t), h(t), p(t) for current t
    fitness[t] <- pars[2] * sum(inputSubset$load * exp(- (t - inputSubset$day)/ pars[3]))
    fatigue[t] <- pars[4] * sum( kh2 * exp(- (t - inputSubset$day) / pars[5]) )
    performance[t] <- pars[1] + fitness[t] - fatigue[t] + initialFitness[t] - initialFatigue[t]
    
  } # Loop index updates (t <- t+1) until t = seriesLength
  
  # Output
  if (returnObject == "performance"){return(performance)}
  if (returnObject == "fitness"){return(fitness)}
  if (returnObject == "fatigue"){return(fatigue)}
  if (returnObject == "all"){
    return(data.frame("day" = 1:seriesLength,
                      "initial_fitness" = initialFitness,
                      "initial_fatigue" = initialFatigue,
                      "fitness" = fitness, "fatigue" = fatigue,
                      "kh2_t" = kh2Dat[, tail(loads$day, 1)],
                      "performance" = performance))
  }
  
} # End function (closing bracket)

# ------------------------------------------------------------------------------
# VDR Model behavior - exploration
# ------------------------------------------------------------------------------

# Constant load values
constantLoad <- data.frame("day" = 0:10,
                           "load" = c(0,rep(1,10)))

# Changing parameter values
parsChanging <- data.frame(c(1, 0.5, 25, 1, 6, 0.2),
                           c(1, 0.5, 25, 1, 6, 0.5),
                           c(1, 0.5, 25, 1, 6, 0.8),
                           c(1, 0.5, 25, 1, 6, 1.2),
                           c(1, 0.5, 25, 1, 6, 1.5),
                           c(1, 0.5, 25, 1, 6, 2))

# Generate predictions
constantPredictions <- lapply(1:6, function(i) simulateVDR(parsChanging[,i], constantLoad))
constantPredictionsFatigue <- sapply(1:6, function(i) cbind(constantPredictions[[i]][,5]))
constantPredictionskh2 <- sapply(1:6, function(i) cbind(constantPredictions[[i]][,6]))

# Generate plots
par(mfrow = c(1,2))
plot(constantPredictionskh2[,6], type = "l", col = "red", ylab = expression(k[h][2](t)),
     xlab = "Day (t)",lty = 2, lwd = 1.5, cex.main = 0.85,
     main = expression("Constant Load, " ~ omega ~ "=" ~ 1 ~ " | " ~ tau[h][2] ~ " varying" ))
lines(constantPredictionskh2[,5], col = "blue", lty = 2, lwd = 1.5)
lines(constantPredictionskh2[,4], col = "green", lty = 2, lwd = 1.5)
lines(constantPredictionskh2[,3], col = "orange", lty = 2, lwd = 1.5)
lines(constantPredictionskh2[,2], col = "pink", lty = 2, lwd = 1.5)
lines(constantPredictionskh2[,1], col = "black", lty = 2, lwd = 1.5)
legend("bottomright", c(expression(tau[h][2] ~ "=" ~ 2), 
                        expression(tau[h][2] ~ "=" ~ 1.5),
                        expression(tau[h][2] ~ "=" ~ 1.2),
                        expression(tau[h][2] ~ "=" ~ 0.8),
                        expression(tau[h][2] ~ "=" ~ 0.5),
                        expression(tau[h][2] ~ "=" ~ 0.2)),
       lty = c(2, 2, 2, 2, 2, 2), lwd = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
       col = c("red", "blue", "green", "orange", "pink", "black"), cex = 0.75)

plot(constantPredictionsFatigue[,6], type = "l", col = "red", ylab = expression("fatigue | " ~ h(t)),
     xlab = "Day (t)",lty = 2, lwd = 1.5, cex.main = 0.6,
     main = expression("Constant Load, " ~ omega ~ "=" ~ 1 ~ " | " ~ tau[h][2] ~ 
                         " varying" ~ " | " ~ tau[h] ~ "=" ~ 6 ~ " | " ~ k[h] ~ "=" ~ 1), cex.main = 0.95)
lines(constantPredictionsFatigue[,5], col = "blue", lty = 2, lwd = 1.5)
lines(constantPredictionsFatigue[,4], col = "green", lty = 2, lwd = 1.5)
lines(constantPredictionsFatigue[,3], col = "orange", lty = 2, lwd = 1.5)
lines(constantPredictionsFatigue[,2], col = "pink", lty = 2, lwd = 1.5)
lines(constantPredictionsFatigue[,1], col = "black", lty = 2, lwd = 1.5)
legend("topleft", c(expression(tau[h][2] ~ "=" ~ 2), 
                    expression(tau[h][2] ~ "=" ~ 1.5),
                    expression(tau[h][2] ~ "=" ~ 1.2),
                    expression(tau[h][2] ~ "=" ~ 0.8),
                    expression(tau[h][2] ~ "=" ~ 0.5),
                    expression(tau[h][2] ~ "=" ~ 0.2)),
       lty = c(2, 2, 2, 2, 2, 2), lwd = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
       col = c("red", "blue", "green", "orange", "pink", "black"), cex = 0.85)


# ------------------------------------------------------------------------------
# 6.8 | Development of synthetic data (via simulation)
# ------------------------------------------------------------------------------

# Import the training load and performance data
# (For purposes of toy example we will create synthetic data)

loads <- data.frame("day" = 0:100, "load" = c(0, rep(c(1, 1.2, 0.5, 1.8, 2, 0.25, 0.7, 0.9, 0, 0.5, 
                                                       1, 0.8, 1.2, 1.3, 0.9, 0, 0, 2, 1.1, 0.5), 5) ) )

# Using the standard model to build some data
mockParameters <- c(100, 1, 22.5, 1.2, 8)       # c(p*, kg, Tg, kh, Th)
mockPerformances <- data.frame("day" = 1:100, "performance" = simulateStandard(mockParameters, loads, returnObject = "performance"))
mockPerformances <- mockPerformances[seq(1, 100, 3),]    # Subset to reduce no. of datapoints

par(mfrow = c(2,1))

par(mar = c(5,5,3,5))
plot(x = mockPerformances$day, y = mockPerformances$performance, type = "l",
     ylab = "Performance [a.u]", main = "Synthetic data", xlab = "Day [a.u]",
     col = "red", lwd = 2)
points(x = mockPerformances$day, y = mockPerformances$performance, col = "red")
par(new = TRUE)
plot(x = loads$day, y = loads$load, type = "h", xlab = "",
     ylab = "", xaxt = "n", yaxt = "n", col = "orange")
mtext("Training load [a.u]", side = 4, line = 3)
axis(side = 4)


# ------------------------------------------------------------------------------
# 6.9 | Fitting the simulated data to the standard and VDR models under NLS via L-BFGS-B
# ------------------------------------------------------------------------------

# Start values for the optimization algorithm
startValsStandard <- c(90, 0.8, 26, 1.5, 11) # c(p*, kg, Tg, kh, Th)
startValsVDR <- c(90, 0.8, 26, 1.5, 11, 0.65) # c(p*, kg, Tg, kh, Th, Th2)
            
# Load required optimisation package (install first if required)
library(optimx)

# Fit the standard model
standardFitted <-  optimx::optimx(par = startValsStandard,      # Starting values
                                  fn = standardObjectiveSS,     # Objective function (RSS)
                                  method = "BFGS",              # Method (see ?optimx)
                                  loads = loads,                # Further arg passed to fn, training loads
                                  perfVals = mockPerformances)  # Further arg passed to fn, performance values

# Fit the VDR model
vdrFitted <-  optimx::optimx(par = startValsVDR,           # Starting values
                             fn = vdrObjectiveSS,          # Objective function (RSS)
                             method = "BFGS",              # Method (see ?optimx)
                             loads = loads,                # Further arg passed to fn, training inputs
                             perfVals = mockPerformances)  # Further arg passed to fn, performance values

# Inspect results
standardFitted
vdrFitted

# ------------------------------------------------------------------------------
# Extra code: Plotting/visualization
# ------------------------------------------------------------------------------

# Compute model predictions
standard_fitted_perf <- simulateStandard(pars = as.numeric(standardFitted[1:5]),
                                                           loads = loads)
vdr_fitted_perf <- simulateVDR(pars = as.numeric(vdrFitted[1:6]),
                               loads = loads)

# Plot results (Figure 6.9)

par(mfrow = c(2,1))
par(mar = c(5,5,3,5))
plot(x = standard_fitted_perf$day, y = standard_fitted_perf$performance, type = "l",
     lwd = 2, ylab = "Performance [a.u]", xlab = "Day [a.u]", main = "Standard Model (Fitted)")
points(x = mockPerformances$day, y = mockPerformances$performance, col = "blue", pch = 16)
par(new = TRUE)
plot(x = standard_fitted_perf$day, y = standard_fitted_perf$fitness, type = "l", xlab = "",
     ylab = "", xaxt = "n", yaxt = "n", col = "green", lty = 2, lwd = 1.5)
lines(x=standard_fitted_perf$day, y = standard_fitted_perf$fatigue, col = "red", lty = 2,
      lwd = 1.5)
mtext("Fitness and Fatigue [a.u]", side = 4, line = 3)
axis(side = 4)
legend("topleft", c("Synthetic data", "Model performance", "Model Fitness", "Model Fatigue"),
       pch = c(16,NA,NA,NA), lty = c(NA,1,2,2), lwd = c(NA,2,1.5,1.5),
       col = c("Blue", "Black", "Green", "Red"), cex = 0.8)

par(mar = c(5,5,3,5))
plot(x = vdr_fitted_perf$day, y = vdr_fitted_perf$performance, type = "l",
     lwd = 2, ylab = "Performance [a.u]", xlab = "Day [a.u]", main = "VDR Model (Fitted)")
points(x = mockPerformances$day, y = mockPerformances$performance, col = "blue", pch = 16)
par(new = TRUE)
plot(x = vdr_fitted_perf$day, y = vdr_fitted_perf$fitness, type = "l", xlab = "",
     ylab = "", xaxt = "n", yaxt = "n", col = "green", lty = 2, lwd = 1.5)
lines(x=vdr_fitted_perf$day, y = vdr_fitted_perf$fatigue, col = "red", lty = 2,
      lwd = 1.5)
mtext("Fitness and Fatigue [a.u]", side = 4, line = 3)
axis(side = 4)
legend("topleft", c("Synthetic data", "Model performance", "Model Fitness", "Model Fatigue"),
       pch = c(16,NA,NA,NA), lty = c(NA,1,2,2), lwd = c(NA,2,1.5,1.5),
       col = c("Blue", "Black", "Green", "Red"), cex = 0.8)