# ******************************************************************************
# Language: R (r-project.org.uk)
# Author: Ben Stephens Hemingway
# License: GNU GPL v3
# Thesis chapter: 6
# Subsection: 6.2.2
# Description: Maximum likelihood estimation (standard, fitness-delay, VDR models)
# R version 4.04 (2021-02-15 "Lost Library Book")
# Package dependencies:
library(optimx)     # General purpose optimization package
# ******************************************************************************

# ------------------------------------------------------------------------------
# 6.12 | Simulation function: Standard model - Sapply approach
# ------------------------------------------------------------------------------

simulateStandard2 <- function(pars, loads, initialPars = c(0,0)){
  
  # Parameters supplied as: pars <- c(p*, k_g, Tau_g, k_h, Tau_h)
  
  # Ancillary function (required)
  # ----------------------------------------------------------------------------
  # Listing 6.10
  convolveTraining <- function(loads, tau){
    
    # Value of t relevant to (eq 6.9) 
    dayt <- length(loads)
    
    # Note that loads[1:dayt] will yield c(w(0), w(1), ... , w(t-1))
    return(sum(loads[1:dayt] * exp(-(dayt:1 / tau))))
  }
  # ----------------------------------------------------------------------------
  
  # Length of the training load series (final day)
  T <- tail(loads$day, 1)
  
  # If initial parameters q_g and q_h are supplied (otherwise evaluates to 0)
  initialFitness <- initialPars[1] * exp(-(1:T) / pars[3])
  initialFatigue <- initialPars[2] * exp(-(1:T) / pars[4])
  
  # Calculate the fitness and fatigue effects (Utilizing sapply function)
  fitness <- pars[2] * 
    base::sapply(1:T, function(t) convolveTraining(loads$load[1:t], 
                                                   pars[3]))
  
  fatigue <- pars[4] * 
    base::sapply(1:T, function(t) convolveTraining(loads$load[1:t], 
                                                   pars[5]))
  
  performance <- pars[1] + initialFitness - initialFatigue + fitness - fatigue
  
  # Return model predicted performance, fitness, and fatigue
  return(data.frame(day = 1:T,
                    performance = performance, 
                    fitness = fitness, 
                    fatigue = fatigue, 
                    load = loads$load[2:(T + 1)]))
}

# ------------------------------------------------------------------------------
# EXTRA: Simulation function: Fitness-delay model - Sapply approach
# ------------------------------------------------------------------------------

simulateFitnessDelay2 <- function(pars, loads, initialPars = c(0,0)){
  
  # Parameters supplied as: pars <- c(p*, k_g, Tau_g, Tau_g2, k_h, Tau_h)
  
  # Ancillary function (required)
  # ----------------------------------------------------------------------------
  convolveTrainingDelay <- function(loads, tau1, tau2){
    
    # Value of t relevant to (eq 6.9) 
    dayt <- length(loads)
    delay <- exp(-(dayt:1 / tau1)) - exp(-(dayt:1 / tau2)) 
    
    # Note that loads[1:dayt] will yield c(w(0), w(1), ... , w(t-1))
    return(sum(loads[1:dayt] * delay))
  }
  
  # Listing 6.10
  convolveTraining <- function(loads, tau){
    
    # Value of t relevant to (eq 6.9) 
    dayt <- length(loads)
    
    # Note that loads[1:dayt] will yield c(w(0), w(1), ... , w(t-1))
    return(sum(loads[1:dayt] * exp(-(dayt:1 / tau))))
  }
  # ----------------------------------------------------------------------------
  
  # Length of the training load series (final day)
  T <- tail(loads$day, 1)
  
  # If initial parameters q_g and q_h are supplied (otherwise evaluates to 0)
  initialFitness <- initialPars[1] * exp(-(1:T) / pars[3])
  initialFatigue <- initialPars[2] * exp(-(1:T) / pars[6])
  
  # Calculate the fitness and fatigue effects (Utilizing sapply function)
  fitness <- pars[2] * 
    base::sapply(1:T, function(t) convolveTrainingDelay(loads$load[1:t], 
                                                        pars[3], pars[4]))
  
  fatigue <- pars[5] * 
    base::sapply(1:T, function(t) convolveTraining(loads$load[1:t], 
                                                   pars[6]))
  
  performance <- pars[1] + initialFitness - initialFatigue + fitness - fatigue
  
  # Return model predicted performance, fitness, and fatigue
  return(data.frame(day = 1:T,
                    performance = performance, 
                    fitness = fitness, 
                    fatigue = fatigue, 
                    load = loads$load[2:(T + 1)]))
}

# ------------------------------------------------------------------------------
# 6.14 | Simulation function: VDR model - Sapply approach
# ------------------------------------------------------------------------------

simulateVDR2 <- function(pars, loads, initialPars = c(0,0)){
  
  # Parameters supplied as: pars <- c(p*, k_g, Tau_g, k_h, Tau_h, Tau_h2)
  
  # Ancillary functions (required)
  # ----------------------------------------------------------------------------
  # Listing 6.10
  convolveTraining <- function(loads, tau){
    
    # Value of t relevant to (eq 6.9) 
    dayt <- length(loads)
    
    # Note that loads[1:dayt] will yield c(w(0), w(1), ... , w(t-1))
    return(sum(loads[1:dayt] * exp(-(dayt:1 / tau))))
  }
  
  # Listing 6.13
  kh2Compute <- function(loads, tau){
    day_i <- length(loads)
    return(sum(loads[1:day_i] * exp(-((day_i-1):0 / tau))))
  }
  # ----------------------------------------------------------------------------
  
  # Length of the training load series supplied = final day (T) for which to compute p\hat(t)
  T <- tail(loads$day, 1)
  
  # If initial parameters q_g and q_h are supplied (otherwise evaluates to 0)
  initialFitness <- initialPars[1] * exp(-(1:T) / pars[3])
  initialFatigue <- initialPars[2] * exp(-(1:T) / pars[5])
  
  # Compute modeled fitness
  fitness <- pars[2] * 
    base::sapply(1:T, function(i) convolveTraining(loads$load[1:i], pars[3]))
  
  # For each i=1:T, compute k_h_2(i) = sum_(j=0)^(i){w(j)e^(-(i-j)/tau_h_2)} 
  # s.t. k_h_2(t) = sum(k_h_2(i)) for i = 1,2,...,(t-1)
  kh2 <- base::sapply(1:T, function(i) kh2Compute(loads$load[1:i], pars[6]))
  
  # Compute modeled fatigue
  fatigue <- pars[4] *
    base::sapply(1:T, function(i) convolveTraining(kh2[1:i], pars[5]))
  
  # Compute modeled performance
  performance <- pars[1] + initialFitness - initialFatigue + fitness - fatigue
  
  # Return model predicted performance, fitness, and fatigue
  return(data.frame(day = 1:T,
                    performance = performance, 
                    fitness = fitness, 
                    fatigue = fatigue,
                    kh2 = kh2,
                    load = loads$load[2:(T + 1)]))
}

# ------------------------------------------------------------------------------
# 6.15 | Objective function (log likelihood): Standard model - Sapply approach
# ------------------------------------------------------------------------------

standardObjectiveLL <- function(pars, loads, perfVals, initial = FALSE, 
                                maximise = FALSE){
  
  # INPUT NOTES:
  # ---------------------------------------------------------------------------
  #         [1] [2] [3] [4] [5] [6]    [7] [8]
  # Pars: c(p*, kg, Tg, kh, Th, sigma)                  initial = FALSE
  # Pars: c(p*, kg, Tg, kh, Th, sigma, qg, qh)          initial = TRUE
  # ---------------------------------------------------------------------------
  
  # Ancillary function (required)
  # ----------------------------------------------------------------------------
  convolveTraining <- function(loads, tau){
    
    # Value of t relevant to (eq 6.9) 
    dayt <- length(loads)
    
    # Note that loads[1:dayt] will yield c(w(0), w(1), ... , w(t-1))
    return(sum(loads[1:dayt] * exp(-(dayt:1 / tau))))
  }
  # ----------------------------------------------------------------------------
  
  finalMeasurement <- tail(perfVals$day, 1)
  
  if (initial == TRUE){
    initFitness <- pars[7] * exp(-(1:finalMeasurement) / pars[3])
    initFatigue <- pars[8] * exp(-(1:finalMeasurement) / pars[5])
  }
  
  # Compute modeled performance from t=1 to t=finalMeasurement
  modFitness <- pars[2] * sapply(1:finalMeasurement, 
                                 function(t) convolveTraining(loads$load[1:t], pars[3]))
  modFatigue <- pars[4] * sapply(1:finalMeasurement, 
                                 function(t) convolveTraining(loads$load[1:t], pars[5]))
  
  if (initial == FALSE){
    modPerformance <- pars[1] + modFitness - modFatigue
  }
  
  if (initial == TRUE){
    modPerformance <- pars[1] + initFitness - initFatigue + modFitness - modFatigue
  }
  
  # Extract modeled performance values on days where measurement exists
  modPerformance <- modPerformance[perfVals$day]
  
  # Compute errors
  errors <- perfVals$performance - modPerformance
  
  if (maximise == FALSE){
    return(-1.0 * sum(dnorm(errors, mean = 0, sd = pars[6], log = TRUE)))
  }
  if (maximise == TRUE){
    return(sum(dnorm(errors, mean = 0, sd = pars[6], log = TRUE)))
  }
}

# ------------------------------------------------------------------------------
# 6.16 | Objective function (log likelihood): VDR model - Sapply approach
# ------------------------------------------------------------------------------

vdrObjectiveLL <- function(pars, loads, perfVals, initial = FALSE, 
                           maximise = FALSE){
  
  # INPUT NOTES:
  # -----------------------------------------------------------------
  #         [1] [2] [3] [4] [5] [6]  [7]    [8] [9]
  # Pars: c(p*, kg, Tg, kh, Th, Th2, sigma)          initial = FALSE
  # Pars: c(p*, kg, Tg, kh, Th, Th2, sigma, qg, qh)  initial = TRUE
  # -----------------------------------------------------------------
  
  # Ancillary functions (required)
  # ----------------------------------------------------------------------------
  convolveTraining <- function(loads, tau){
    
    # Value of t relevant to (eq 6.9) 
    dayt <- length(loads)
    
    # Note that loads[1:dayt] will yield c(w(0), w(1), ... , w(t-1))
    return(sum(loads[1:dayt] * exp(-(dayt:1 / tau))))
  }
  
  kh2Compute <- function(loads, tau){
    day_i <- length(loads)
    return(sum(loads[1:day_i] * exp(-((day_i-1):0 / tau))))
  }
  # ----------------------------------------------------------------------------
  
  finalMeasurement <- tail(perfVals$day, 1)
  
  if (initial == TRUE){
    initFitness <- pars[8] * exp(-(1:finalMeasurement) / pars[3])
    initFatigue <- pars[9] * exp(-(1:finalMeasurement) / pars[5])
  }
  
  # Compute modeled performance from t=1 to t=finalMeasurement
  fitness <- pars[2] * sapply(1:finalMeasurement, 
                                 function(t) convolveTraining(loads$load[1:t], pars[3]))
  kh2 <- sapply(1:finalMeasurement, function(i) kh2Compute(loads$load[1:i], pars[6]))
  fatigue <- pars[4] * sapply(1:finalMeasurement, 
                                 function(t) convolveTraining(kh2[1:t], pars[5]))
  
  if (initial == FALSE){
    performance <- pars[1] + fitness - fatigue
  }
  
  if (initial == TRUE){
    performance <- pars[1] + initFitness - initFatigue + fitness - fatigue
  }
  
  # Extract modeled performance values on days where measurement exists
  performance <- performance[perfVals$day]
  
  # Compute errors
  errors <- perfVals$performance - performance
  
  if (maximise == FALSE){
    return(-1.0 * sum(dnorm(errors, mean = 0, sd = pars[7], log = TRUE)))
  }
  if (maximise == TRUE){
    return(sum(dnorm(errors, mean = 0, sd = pars[7], log = TRUE)))
  }
}

# ------------------------------------------------------------------------------
# EXTRA | Objective function (log likelihood): Fitness-delay model - Sapply approach
# ------------------------------------------------------------------------------

fitnessDelayObjectiveLL <- function(pars, loads, perfVals, initial = FALSE, 
                                    maximise = FALSE){
  
  # INPUT NOTES:
  # ---------------------------------------------------------
  #         [1] [2] [3] [4]  [5] [6] [7]    [8] [9]
  # Pars: c(p*, kg, Tg, Tg2, kh, Th, sigma)          initial = FALSE
  # Pars: c(p*, kg, Tg, Tg2, kh, Th, sigma, qg, qh)  initial = TRUE
  # ----------------------------------------------------------
  
  # Ancillary function (required)
  # ----------------------------------------------------------------------------
  convolveTrainingDelay <- function(loads, tau1, tau2){
    
    # Value of t relevant to (eq 6.9) 
    dayt <- length(loads)
    delay <- exp(-(dayt:1 / tau1)) - exp(-(dayt:1 / tau2)) 
    
    # Note that loads[1:dayt] will yield c(w(0), w(1), ... , w(t-1))
    return(sum(loads[1:dayt] * delay))
  }
  
  convolveTraining <- function(loads, tau){
    
    # Value of t relevant to (eq 6.9) 
    dayt <- length(loads)
    
    # Note that loads[1:dayt] will yield c(w(0), w(1), ... , w(t-1))
    return(sum(loads[1:dayt] * exp(-(dayt:1 / tau))))
  }
  # ----------------------------------------------------------------------------
  
  finalMeasurement <- tail(perfVals$day, 1)
  
  if (initial == TRUE){
    initFitness <- pars[8] * exp(-(1:finalMeasurement) / pars[3])
    initFatigue <- pars[9] * exp(-(1:finalMeasurement) / pars[6])
  }
  
  # Compute modeled performance from t=1 to t=finalMeasurement
  modFitness <- pars[2] * 
    sapply(1:finalMeasurement, function(t) convolveTrainingDelay(loads$load[1:t], pars[3], pars[4]))
  modFatigue <- pars[5] * 
    sapply(1:finalMeasurement, function(t) convolveTraining(loads$load[1:t], pars[6]))
  
  if (initial == FALSE){
    modPerformance <- pars[1] + modFitness - modFatigue
  }
  
  if (initial == TRUE){
    modPerformance <- pars[1] + initFitness - initFatigue + modFitness - modFatigue
  }
  
  # Extract modeled performance values on days where measurement exists
  modPerformance <- modPerformance[perfVals$day]
  
  # Compute errors
  errors <- perfVals$performance - modPerformance
  
  if (maximise == FALSE){
    return(-1.0 * sum(dnorm(errors, mean = 0, sd = pars[7], log = TRUE)))
  }
  if (maximise == TRUE){
    return(sum(dnorm(errors, mean = 0, sd = pars[7], log = TRUE)))
  }
  
}

# ------------------------------------------------------------------------------
# 6.17 | Objective function (log likelihood): VDR model - Sapply approach
# ------------------------------------------------------------------------------

loads <- data.frame("day" = 0:100, "load" = c(0, rep(c(1, 1.2, 0.5, 1.8, 2, 0.25, 0.7, 0.9, 0, 0.5, 
                                                       1, 0.8, 1.2, 1.3, 0.9, 0, 0, 2, 1.1, 0.5), 5) ) )

# Using the standard model to build the synthetic performance data
mockParameters <- c(100, 1, 22.5, 1.2, 8, 0.5)       # c(p*, kg, Tg, kh, Th, Th2)
mockPerformances <- simulateVDR2(mockParameters, loads)[,1:2]
mockPerformances_noerror <- mockPerformances
# Add some random error from a Gaussian distribution
sigma <- 1
set.seed(1)
mockPerformances$performance <- mockPerformances$performance + rnorm(100, 0, sigma)
# Subset to reduce number of data points
mockPerformances <- mockPerformances[seq(1, 100, 3),]

# ------------------------------------------------------------------------------
# Plot (figure 6.10)
# ------------------------------------------------------------------------------
par(mfrow = c(1,2), mai = c(0.9, 0.8, 0.8, 0.8))
plot(x=mockPerformances$day, y= mockPerformances_noerror$performance[seq(1,100,3)], col = "blue", pch = 1,
     ylab = "Performance (a.u)", xlab = "Day", main = expression("Without Gaussian error"), cex.main = 0.75)
legend("bottomright", "Simulated data (VDR model)", pch = 1, col = "blue", cex = 0.75)

plot(x=mockPerformances$day, y=mockPerformances$performance,
     ylab = "Performance (a.u)", 
     xlab = "Day",
     pch = 1, col = "red", main = expression("+ random Gaussian error"~(sigma~"="~1)), cex.main = 0.75)
legend("bottomright", "Simulated data (VDR model)", pch = 1, col = "red", cex = 0.75)


# ------------------------------------------------------------------------------
# 6.18 | Fitting the standard and VDR models to the simulated data via MLE & L-BFGS-B
# ------------------------------------------------------------------------------

# Fit the standard model
startPars <- list("standard" = c(95, 0.85, 26, 1.4, 5, 0.6),    #(p*, kg, Tg, kh, Th, sigma)
                   "vdr" = c(95, 0.85, 26, 1.4, 5, 1, 1.3))     #(p*, kg, Tg, kh, Th, Th2, sigma)

standardFitted <- optimx(par = startPars[["standard"]], 
                         fn = standardObjectiveLL,
                         lower = c(60, 0.1, 1, 0.1, 1, 0.1), 
                         upper = c(200, 3, 50, 3, 50, 10),
                         method = "L-BFGS-B",
                         loads = loads,
                         perfVals = mockPerformances,
                         control = list(maxit = 10000))

# Fit the VDR model
vdrFitted <- optimx(par = startPars[["vdr"]], 
                    fn = vdrObjectiveLL,
                    lower = c(60, 0.1, 1, 0.1, 1, 0.1, 0.1), 
                    upper = c(200, 3, 50, 3, 50, 5, 10),
                    method = "L-BFGS-B",
                    loads = loads,
                    perfVals = mockPerformances,
                    control = list(maxit = 10000))

# Compute the associated predictions

standardPredicted <- simulateStandard2(pars = as.numeric(standardFitted[1:5]), loads)
vdrPredicted <- simulateVDR2(pars = as.numeric(vdrFitted[1:6]), loads)

# Plot the results of the fitting process

plot(mockPerformances, ylab = "Performance [a.u]", xlab = "Day", 
     main = "Fitted Model's (Standard FFM and VDR)", lwd = 2)
lines(standardPredicted$performance, col = "blue", lty = 2, lwd = 1.5)
lines(vdrPredicted$performance, col = "red", lwd = 1, lty = 1)
legend("topleft", c("Standard Model (fitted)", "VDR Model (fitted)", "Measured Performance"), 
       lty = c(1,2,NA), pch = c(NA,NA,1), lwd = c(1,2,1.5),
       col = c("red", "blue", "black"), cex = 0.85)