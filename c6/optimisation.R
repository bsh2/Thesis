# ******************************************************************************
# Language: R (r-project.org.uk)
# Author: Ben Stephens Hemingway
# License: GNU GPL v3
# Thesis chapter: 6
# Subsection: 6.2.4
# Description: Model fitting (optimisation)
# R version 4.04 (2021-02-15 "Lost Library Book")
# Package dependencies:
library(optimx)
library(GA)
library(pso)
library(cmaes)
library(DEoptim)
# ******************************************************************************

# ------------------------------------------------------------------------------
# IMPORT REQUIRED MODELLING FUNCTIONS
# ------------------------------------------------------------------------------
# Data and functions required
source("maximum_likelihood.R")
# Clear plots and remove unnecessary objects from sourcing the file
# (Leaves required performance data and loads only)
dev.off()
rm(mockPerformances, mockParameters, sigma, standardFitted, vdrFitted,
   simulateVDR2, simulateStandard2, simulateFitnessDelay2, fitnessDelayObjectiveLL,
   vdrPredicted, standardPredicted, startPars)


# ------------------------------------------------------------------------------
# 6.26 | Fit the standard model, via multistart function (optimx)
#        MLE (L-BFGS-B)
# ------------------------------------------------------------------------------

# Develop matrix of parameter values (you could do this randomly)
pars <- as.matrix(expand.grid("p0" = seq(80, 120, length.out = 2), "kg" = seq(0.5, 2, length.out = 2),
              "Tg" = seq(5, 40, length.out = 2), "kh" = seq(0.5, 2, length.out = 2),
              "Th" = seq(5, 40, length.out = 2), "sigma" = seq(0.1, 1, length.out = 2)))

# Apply the multistart algorithm under the "L-BFGS-B" algorithm in optimx
standard_model1 <- optimx::multistart(parmat = pars,
                                     fn = standardObjectiveLL,
                                     method = "L-BFGS-B",
                                     # ORDER: c(p0, kg, Tg, kh, Th, sigma)
                                     lower = c(60, 0.1, 1, 0.1, 1, 0.1), 
                                     upper = c(200, 3, 50, 3, 50, 10),
                                     loads = loads,
                                     perfVals = mockPerformances_noerror)
head(standard_model)

# ------------------------------------------------------------------------------
# 6.27 | Fit the standard model, MLE (DEOptim)
# ------------------------------------------------------------------------------

# DEOPTIM PACKAGE
standard_model2 <- DEoptim::DEoptim(standardObjectiveLL,
                                   lower = c(60, 0.1, 1, 0.1, 1, 0.1),
                                   upper = c(200, 3, 50, 3, 50, 10),
                                   DEoptim.control(
                                     iter = 1000,
                                   ),
                                   loads = loads,
                                   perfVals = mockPerformances_noerror)

# ------------------------------------------------------------------------------
# 6.28 | Fit the standard model, MLE (GA, CMA-ES, Particle Swarm)
# ------------------------------------------------------------------------------

# GA
standard_model3 <- GA::ga(type = "real-valued",
                         fitness = standardObjectiveLL,
                         perfVals = mockPerformances_noerror,
                         loads = loads,
                         lower = c(60, 0.1, 1, 0.1, 1, 0.1),
                         upper = c(200, 3, 50, 3, 50, 10),
                         maxiter = 1000,
                         monitor = TRUE,
                         optim = TRUE,
                         maximise = TRUE)

# PSO
standard_model4 <- pso::psoptim(par = c(NA,NA,NA,NA,NA,NA,NA),
                               fn = standardObjectiveLL,
                               lower = c(60, 0.1, 1, 0.1, 1, 0.1),
                               upper = c(200, 3, 50, 3, 50, 10),
                               control = list(maxit = 1000),
                               loads = loads,
                               perfVals = mockPerformances_noerror)

# CMAES
standard_model5 <- cmaes::cma_es(par = c(95, 1, 28, 2, 5, 0.1), # guess
                                fn = standardObjectiveLL,
                                lower = c(60, 0.1, 1, 0.1, 1, 0.1),
                                upper = c(200, 3, 50, 3, 50, 10),
                                control = list(maxit = 1000),
                                loads = loads,
                                perfVals = mockPerformances_noerror)
