# ******************************************************************************
# Language: R (r-project.org.uk)
# Author: Ben Stephens Hemingway
# License: GNU GPL v3
# Thesis chapter: 6
# Subsection: 6.1.1.1
# Description: Solve and fit original FFM and nonlinear ODE systems
# R version 4.04 (2021-02-15 "Lost Library Book")
# Package dependencies:
library(deSolve)  # ODE solver
library(optimx)   # First and second-order optimisation algorithms
library(GA)       # Genetic algorithm (global optimiser)
# ******************************************************************************

# ------------------------------------------------------------------------------
# 6.19 | Defining the standard model ODE system as an R function
# ------------------------------------------------------------------------------

banisterSystem <- function(t, y, parms, currentLoad){
  r = c()
  r[1] <- (parms[2]*currentLoad) - ((1/parms[3]) * y["G"])
  r[2] <- (parms[4]*currentLoad) - ((1/parms[5]) * y["H"])
  return(list(r))
}

# ------------------------------------------------------------------------------
# 6.20| Applying a numerical integrator to develop a model simulation function 
#       from the original system of ODEs (standard model)
# ------------------------------------------------------------------------------

banisterSimulate <- function(parms, loads){
  
  dat <- loads
  dat$G <- c(rep(0, length(loads$day)))
  dat$H <- c(rep(0, length(loads$day)))
  dat$pHat <- c(rep(0, length(loads$day)))
  
  # Solve model numerically at each time point (j)
  for (j in 1:length(dat$day)){
    currentLoad <- dat$load[j]
    if (j < 3){
      stateInit <- c(G = parms[6],     # Fitness (initial condition)
                     H = parms[7])     # Fatigue (initial condition)
    } else{
      # Initialize based on previous value (j-1)
      stateInit <- c(G = dat$G[j-1],    # Fitness on previous day 
                     H = dat$H[j-1])    # Fatigue on previous day
    }
    
    # Vector of time-steps up to next day (j+1)
    t <- 0:1
    
    # Solve model for current time point (j)
    out = ode(y = stateInit, times = t, func = banisterSystem, parms = parms, 
              method = c("lsoda"), currentLoad = currentLoad)
    
        # Other integrator methods available include lsoda, euler, 
    
    # Extract solutions (j)
    if (j < 3){
      dat$G[j] <- unname(out[1,2])
      dat$H[j] <- unname(out[1,3])
    } else{
      dat$G[j] <- unname(out[2,2])
      dat$H[j] <- unname(out[2,3])
    }
    
    # Modelled performance at time j
    dat$pHat[j] <- parms[1] + dat$G[j] - dat$H[j]
    
  } # Update loop (j+1) and repeat
  
  return(dat)
}

# ------------------------------------------------------------------------------
# 6.21 | RSS wrapper to facilitate NLS estimation of FFM system parameters & ICS
# ------------------------------------------------------------------------------

banisterObjective <- function(parms, perfVals, loads){
  
  dat <- banisterSimulate(parms, loads)
  datSubset <- dat[c((perfVals$day) + 1), ]   # +1 due to day zero first row
  
  RSS <- sum((datSubset$pHat - perfVals$performance)^2)
  return(RSS)
}

# ------------------------------------------------------------------------------
# 6.22 | Solving the original ODE system by numerical methods and fitting to
#        data (from listing 6.17), via NLS solver (L-BFGS-B)
# ------------------------------------------------------------------------------

# Data required (place maximum_likelihood.R in source file location)
source("maximum_likelihood.R")
# Clear plots and remove unnecessary objects from sourcing the file
# (Leaves required performance data and loads only)
dev.off()
rm(mockPerformances_noerror, standardFitted, standardPredicted, startPars,
   vdrFitted, vdrPredicted, mockParameters, sigma,
   simulateStandard2, simulateFitnessDelay2, simulateVDR2,
   fitnessDelayObjectiveLL, standardObjectiveLL, vdrObjectiveLL)


# Fitting example: listing 6.22
  # Create some starting values for the search

  startParameters <- c(100, 0.9, 26, 1.2, 5, 2, 2)
  
  fittedModel <- optimx(par = startParameters,
                 fn = banisterObjective,
                 method = "L-BFGS-B",
                 lower = c(50, 0.1, 1, 0.1, 1, 0, 0),
                 upper = c(150, 3, 50, 3, 50, 20, 20),
                 perfVals = mockPerformances,  # Same performance data as listing 6.17
                 loads = loads)                # Same loads data as listing 6.17

# ------------------------------------------------------------------------------
# TURNER MODEL (FUNCTIONS) - MODIFICATIONS TO ABOVE CODE (Inc. Listing 23)
# ------------------------------------------------------------------------------

  # System function
  turnerSystem <- function(t, y, parms, currentLoad){
    r = c()
    r[1] <- (parms[2]*currentLoad) - ((1/parms[3]) * y["G"]^parms[8]) # Modification 1
    r[2] <- (parms[4]*currentLoad) - ((1/parms[5]) * y["H"]^parms[9]) # Modification 2
    return(list(r))
  }
  
  # Simulation function
  turnerSimulate <- function(parms, loads){
    
    dat <- loads
    dat$G <- c(rep(0, length(loads$day)))
    dat$H <- c(rep(0, length(loads$day)))
    dat$pHat <- c(rep(0, length(loads$day)))
    
    # Solve model numerically at each time point (j)
    for (j in 1:length(dat$day)){
      currentLoad <- dat$load[j]
      if (j < 3){
        stateInit <- c(G = parms[6],     # Fitness (initial condition)
                       H = parms[7])     # Fatigue (initial condition)
      } else{
        # Initialize based on previous value (j-1)
        stateInit <- c(G = dat$G[j-1],    # Fitness on previous day 
                       H = dat$H[j-1])    # Fatigue on previous day
      }
      
      # Vector of time-steps up to next day (j+1)
      t <- 0:1
      
      # Solve model for current time point (j)
      out = ode(y = stateInit, times = t, func = turnerSystem, parms = parms, 
                method = c("lsoda"), currentLoad = currentLoad)
      
      # Other integrator methods available include lsoda, euler, 
      
      # Extract solutions (j)
      if (j < 3){
        dat$G[j] <- unname(out[1,2])
        dat$H[j] <- unname(out[1,3])
      } else{
        dat$G[j] <- unname(out[2,2])
        dat$H[j] <- unname(out[2,3])
      }
      
      # Modelled performance at time j
      dat$pHat[j] <- parms[1] + dat$G[j] - dat$H[j]
      
    } # Update loop (j+1) and repeat
    
    return(dat)
  }
  
  # Objective function
  turnerObjective <- function(parms, perfVals, loads){
    
    dat <- turnerSimulate(parms, loads)
    datSubset <- dat[c((perfVals$day)+1), ]   # +1 due to day zero first row
    
    RSS <- sum((datSubset$pHat - perfVals$performance)^2)
    return(RSS)
  }
  
  # Fitting example: 
  # Create some starting values for the search
  startParameters <- c(100, 0.9, 26, 1.2, 5, 0, 0, 1, 1) # Modification 3
  
  # Need to bound Turner as weird things happen when alpha or beta nears <=0
  constraints <- data.frame(lower = c(50, 0.1, 1, 0.1, 1, 0.1, 0.1, 0.4, 0.4),
                            upper = c(200, 10, 50, 10, 50, 20, 20, 3, 3))
  
  fittedModel <- optimx(par = startParameters,
                        fn = turnerObjective,
                        method = "L-BFGS-B",
                        lower = constraints$lower,
                        upper = constraints$upper,
                        perfVals = mockPerformances,      # Same performance data as listing 6.17
                        loads = loads)                    # Same loads data as listing 6.17
  

# ------------------------------------------------------------------------------
# TURNER MODEL (PLOTS) - Constant Load
# ------------------------------------------------------------------------------
  
  # The training load that causes optimal (max) performance under constant load
  w_optimal <- function(alpha, beta, kg, Tau_g, kh, Tau_h){
    out <- (((alpha/beta)^(alpha * beta)) *
           (((kh * Tau_h)^(alpha)) / ((kg * Tau_g)^(beta))))^(1/(beta-alpha))
    return(out)
  }
  
  # The training load value that would result in no performance improvement
  w_max <- function(kg, Tau_g, kh, Tau_h, alpha, beta){
    return((((kh*Tau_h)^(alpha)) / 
              ((kg*Tau_g)^(beta)))^(1/(beta-alpha)))
  }
  
  # Steady state performance under constant load
  p_tilde <- function(p_0, kg, Tau_g, kh, Tau_h, alpha, beta, w){
    return(p_0 + (kg*Tau_g*w)^(1/alpha) - (kh*Tau_h*w)^(1/beta))
  }
  
  performance <- p_tilde(p_0 = 155, kg = 0.10, Tau_g = 61, kh = 0.12, Tau_h = 5.5, 
                         alpha = 1.16, beta = 0.85, w = seq(0,750,1))
  max_performance <- max(performance)
  w_o <- w_optimal(alpha = 1.16, beta = 0.85, kg = 0.10, Tau_g = 61, kh = 0.12, Tau_h = 5.5)
  w_m <- w_max(kg = 0.10, Tau_g = 61, kh = 0.12, Tau_h = 5.5, alpha = 1.16, beta = 0.85)

  # Figure 6.12
  plot(x = seq(0,750, 1), y = performance, lwd = 2,
       type = "l", col = "darkblue", ylab = expression(tilde(p)~"[performance]"),
       xlab = expression(omega ~ "[Constant load]"))
  abline(h = 155, lty = 2) # Line at p* (Baseline)
  abline(v = w_o, lty = 3, col = "forestgreen") # Line at w = w_optimal
  abline(v = w_m, lty = 4, col = "red") # Line at w = w_max
  points(x = w_o, y = max_performance, col = "orange", pch = 16) # Maximum performance
  legend("topright", c(expression(tilde(p[omega])), expression(p~"*"),
                       expression(omega["optimal"]), expression(omega["max"]),
                       expression("performance"["max"])),
         lty = c(1, 2, 3, 4, NA), pch = c(NA, NA, NA, NA, 16),
         col = c("darkblue", "black", "forestgreen", "red", "orange"), cex = 0.9)
  
  # Figure 6.13
  w <- seq(0,750,1); kg = 0.1; kh = 0.12; Tau_g = 61; Tau_h = 5.5; alpha = 1.16; beta = 0.85
  fitness <- (kg*Tau_g*w)^(1/alpha)
  fatigue <- (kh*Tau_h*w)^(1/beta)
  
  plot(w, fitness, type = "l", col = "blue", xlab = expression(omega ~ "[Constant load]"),
       ylab = expression(tilde(g)~"and"~tilde(h)~"[a.u]"))
  lines(w, fatigue, col = "red", lty = 2)
  points(w_m, fatigue[w_m]) # Point at which the lines cross and fatigue overwhelm's fitness
  abline(v = w_o, lty = 3, col = "forestgreen") # Line at w_optimal
  abline(v = w_m, lty = 4, col = "black") # Line at w = w_max
  legend("topleft", c(expression(tilde(g)), expression(tilde(h)), expression(omega[optimal]),
                       expression(omega[max])),
         lty = c(1,2,3,4), col = c("blue", "red", "forestgreen", "black"))

# ------------------------------------------------------------------------------
# 6.25 | Fitting Turner's model to the synthetic data via genetic algorithm
# ------------------------------------------------------------------------------

# Packages required
require(parallel)
require(doParallel)
require(GA)
  
# Need to bound Turner as weird things happen when alpha or beta nears <=0
constraints <- data.frame(lower = c(50, 0.1, 1, 0.1, 1, 0.1, 0.1, 0.4, 0.4),
                          upper = c(200, 10, 50, 10, 50, 20, 20, 3, 3))

# You must first adjust the objective function above!!!
# To change a max problem to a min problem, just multiply the objective function by −1

# Objective function: listing 6.20
turnerObjective <- function(parms, perfVals, loads){
  
  dat <- turnerPredict(parms, loads)
  datSubset <- dat[c((perfVals$day)+1), ]   # +1 due to day zero first row
  
  RSS <- sum((datSubset$pHat - perfVals$performance)^2)
  return(-RSS)
  # GA is a maximiser so we will look to maximise the negative RSS
}

fittedModel <- GA::ga(type = "real-valued",      # Type
                      fitness = turnerObjective, # Objective function
                      perfVals = mockPerformances, # Measured performance values
                      loads = loads,               # Training loads
                      lower = constraints$lower, # Lower bound
                      upper = constraints$upper, # Upper bound
                      maxiter = 1000,     # Maximum number of iterations
                      monitor = TRUE,  # Console output during process
                      popSize = 90,       # Population size
                      optim = FALSE,       # Include local search (Y/N)
                      # if (optim = TRUE) the following reflect the local search args
                      # optimArgs = list(method = c("L-BFGS-B"),
                      #                  poptim = 0.1, # [0,1] probability of performing a local search
                      #                  pressel = 0.5, # [0,1] pressure selection (default 0.5)
                      #                  lower = constraints$lower,
                      #                  upper = constraints$upper,
                      #                 control = list(maxit = 1500)),
                      elitism = 7.5,                    # Survival at each generation (top %)
                      selection = gareal_tourSelection, # Tournament selection
                      crossover = gareal_blxCrossover,  # BLX-alpha crossover
                      mutation = gareal_rsMutation,     # Random Gaussian mutation
                      run = 150,                        # Halt value
                      parallel = "multicore",           # Multi-platform parallelisation
                      seed = 12345                      # Seed for replication later
                ) # End of GA call

  