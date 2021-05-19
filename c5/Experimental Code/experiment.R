# Step 0: Import functions and setup environment

  # Required functions

    source("functions/standardModel.R") # Loss, optim, compute functions
    source("functions/calvertModel.R")  # Loss, optim, compute functions
    source("functions/metrics.R")       # R squared, RMSE, MAPE functions (fit)
    source("functions/gridSearchOptim.R")     # Grid Search + Local (optim)
    
  # Required Packages
    
    library(stats)
    library(parallel)
    library(snow)
    library(doSNOW)
    
  # Grid Size
    
    gridDimStandard <- 10  # Number of total sets = gridDimStandard^(nPars = 5)
    gridDimCalvert <- 7   # Number of total sets = gridDimCalvert^(nPars = 6)
    
    gridDivStandard <- 2
    gridDivCalvert <- 7
    
    # Note:
    #     Due to memory limitations, the grid needs to be split and evaluated
    #     in sections. Hence we must perform the check below
    
    if (gridDimStandard^5 %% gridDivStandard != 0){
      stop("Non-divisible grid")
    }
    
    if (gridDimCalvert^6 %% gridDivCalvert != 0){
      stop("Non-divisible grid")
    }
    
    dir.create("output")
    dir.create("output/standard")
    dir.create("output/calvert")
    

for (freqIndex in 1:3){
  
  # Set up a training load series
  
    loads <- read.csv("loads.csv")
    colnames(loads) <- c("load")
    loads <- as.numeric(loads$load)
  
  # Establish a deterministic (true) state for the standard model
  
    # Set the deterministic true parameters and compute the reference state
    standardPars <- c(100, 0.72, 28.5, 1.2, 8.6)
    standardFFM <- standardCompute(standardPars, loads)
  
  # Measurement frequency
  
    criterionFrequency <- freqIndex
    # Note to self: As measurement frequency decreases, there should still be an
    # RSS = 0 state, if there is no additional noise, so it will be interesting
    # to note if this gets harder to find.
  
  # Organize the true state data ready for fitting (Subset by measurement freq)
  
  standardIndex <- seq(1, length(standardFFM$performance), 
                       by = criterionFrequency)
  standardCriterion <- standardFFM$performance[standardIndex]
  
# Establish a deterministic (true) state for the calvert (delay) model
  
  # Set the deterministic true parameters and compute the reference state
  calvertPars <- c(100, 0.72, 32.5, 4.3, 1.05, 8.6)
  calvertFFM <- calvertCompute(calvertPars, loads)
  
  
  # Organize the true state data ready for fitting (Subset by measurement freq)
  calvertIndex <- seq(1, length(calvertFFM$performance), 
                  by = criterionFrequency)
  calvertCriterion <- calvertFFM$performance[calvertIndex]
  
# Step 4: Set up optimization conditions
  
  # c(p*, kg, Tg, kh, Th)
  standardConstraints <- data.frame("lower" = c(60, 0.01, 1, 0.01, 1),
                                    "upper" = c(140, 5, 50, 5, 50))
  
  # c(p*, kg, Tg1, Tg2, kh, Th)
  calvertConstraints <- data.frame("lower" = c(60, 0.01, 1, 1, 0.01, 1),
                                   "upper" = c(140, 5, 50, 50, 5, 50))
  
# Step 5: Set up the grid(s) of starting values and run optimization
  
  dir.create(paste0("output/standard/", freqIndex))
  dir.create(paste0("output/calvert/", freqIndex))

  # Standard Model
  
    # Grid Generate, Evaluate Starting Loss, and call Optim() from the points
    standardGrid <- gridStandard(freqIndex)
    rm(standardGrid)
    gc()
    
  # Calvert Model
    
    # Grid Generate, Evaluate Starting Loss, and call Optim() from the points
    calvertGrid <- gridCalvert(freqIndex)
    rm(calvertGrid)
    gc()

} # End Loop