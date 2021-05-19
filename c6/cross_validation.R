# ******************************************************************************
# Language: R (r-project.org.uk)
# Author: Ben Stephens Hemingway
# License: GNU GPL v3
# Thesis chapter: 6
# Subsection: 6.2.5
# Description: Cross-validation (expanding-window approach)
# R version 4.04 (2021-02-15 "Lost Library Book")
# Package dependencies: 
library(optimx)
library(caret)      # Use createTimeSlices function to derive our splits
library(RcppAlgos)  # Use comboSample function to derive starting points
library(foreach)    # For training and validating splits in parallel
library(parallel)   # Parallel backend
library(doSNOW)     # Parallel backend
# ******************************************************************************

# ------------------------------------------------------------------------------
# Data and functions required
# ------------------------------------------------------------------------------
source("maximum_likelihood.R")
# Clear plots and remove unnecessary objects from sourcing the file
# (Leaves required performance data and loads only)
dev.off()
rm(mockPerformances, mockParameters, sigma, standardFitted, vdrFitted,
   simulateStandard2, simulateFitnessDelay2, fitnessDelayObjectiveLL,
   vdrPredicted, standardPredicted, startPars, standardObjectiveLL)

# ------------------------------------------------------------------------------
# 6.29 | Developing mock data used to demonstrate the CV implementation
# ------------------------------------------------------------------------------

block1_loads <- rep(c(1, 1.2, 0.5, 1.8, 2, 0.25, 0.7,0.9, 0, 0.5, 1, 0.8, 
                      1.2, 1.3, 0.9, 0, 0, 2, 1.1, 0.5), 5)

block2_loads <- round(abs(rnorm(50, 1.3, 0.5)),1)

loads <- data.frame("day" = 0:150,
                    "load" = c(0, block1_loads, block2_loads),
                    "block" = c(rep(1, 101), rep(2, 50)))

rm(block1_loads, block2_loads)

# Generating two blocks of associated performance data
true_parameters <- c(95, 0.85, 26, 1.2, 5, 1) #c(p*,kg,Tg,kh,Th,Th2)
performances <- simulateVDR2(true_parameters, loads = loads)
performances$block <- loads$block[2:151]

# Reduce measurement frequency and regularity
set.seed(101)
performances <- performances[c(1, sort(as.integer(sample(2:149, 56), replace = FALSE), 
                                       decreasing = FALSE), 150), ]
  
# Plot mock data
par(mfrow = c(2,1))
plot(x = performances$day, y = performances$performance, type = "n",  pch = 1,
     col = "blue", ylab = "Performance [a.u]", xlab = "Day",
     main = "Synthetic (simulated) performance", cex.main = 0.95, cex = 0.95,
     xlim = c(0,150))
points(x = performances[performances$day <= 100,"day"], 
       y = performances[performances$day <= 100,"performance"], col = "blue")
points(x = performances[performances$day > 100,"day"], 
       y = performances[performances$day > 100,"performance"], col = "brown2")
abline(v = 100, lty = 2, lwd = 1.5, col = "green")
text(x = 128, y = 95, labels = c("Block 2"), col = "brown2")
text(x = 50, y = 95, labels = c("Block 1"), col = "blue")
plot(x = loads$day, y = loads$load, type = "n",
     ylab = "Training load [a.u]", xlab = "Day", cex.main = 0.95,
     main = "Synthetic training loads", ylim = c(0, 3.5), col = "blue")
lines(x=loads[loads$day <= 100, "day"], y= loads[loads$day <= 100, "load"], type = "h", col = "blue")
lines(x=loads[loads$day > 100, "day"], y= loads[loads$day > 100, "load"], type = "h", col = "brown2")
abline(v = 100, lty = 2, lwd = 1.5, col = "green")
text(x=50, y= 3, labels = c("Block 1"), col = "blue")
text(x=128, y = 3, labels = c("Block 2"), col = "brown2")
  
# Develop some bounds for the parameter space
bounds <- data.frame("lower" = c(80, 0.1, 2, 0.1, 2, 0.2, 0.2),
                       "upper" = c(120, 3, 50, 3, 50, 3, 5))
  
# Collate the data (input format for the CV function)
dat <- data.frame("day" = loads$day,
                  "load" = loads$load,
                  "performance" = rep(NA, length(loads$load)),
                  "block" = loads$block)
dat$performance[performances$day + 1] <- performances$performance
rm(loads, performances)

# ------------------------------------------------------------------------------
# 6.30 | MAPE function
# ------------------------------------------------------------------------------

mape <- function(measured, predicted){
  return(mean(abs((measured - predicted)/measured))*100)
}

# ------------------------------------------------------------------------------
# 6.31 | Function to develop expanding window splits
# ------------------------------------------------------------------------------

generate_splits <- function(initialWindow, testHorizon, expandRate, dat){
  # Convert supplied arguments for CV into 'days'
  initialWindow = round(length(dat$day) * initialWindow/100, 0)
  testHorizon = round(length(dat$day) * testHorizon/100, 0)
  expandRate = round(length(dat$day) * expandRate/100, 0)
  # Create splits
  splits <- createTimeSlices(dat$day,
                             initialWindow = initialWindow,
                             horizon = testHorizon,
                             fixedWindow = FALSE,
                             skip = expandRate)
  return(splits)
}

# ------------------------------------------------------------------------------
# 6.32 | Developing a function to create a random grid of starting parameters
#        across the bounds [l,u] to be used for iterative model fitting from
#        multiple starting points
# ------------------------------------------------------------------------------


create_grid <- function(bounds, nStarts, initial){
  
  if (initial == FALSE){
    set.seed(101)
    parmat <- data.frame(
      "p0" = comboSample(seq(bounds$lower[1] + 1, bounds$upper[1] - 1, 2.5), 
                         m = 1, n = nStarts),
      "kg" = comboSample(seq(bounds$lower[2] + 0.1, bounds$upper[2] - 0.1, 0.1), 
                         m = 1, n = nStarts),
      "Tg" = comboSample(seq(bounds$lower[3] + 1, bounds$upper[3] - 1, 1.5), 
                         m = 1, n = nStarts),
      "kh" = comboSample(seq(bounds$lower[4] + 0.1, bounds$upper[4] - 0.1, 0.1), 
                         m = 1, n = nStarts),
      "Th" = comboSample(seq(bounds$lower[5] + 1, bounds$upper[5] - 1, 1.5), 
                         m = 1, n = nStarts),
      "Th2" = comboSample(seq(bounds$lower[6] + 0.1, bounds$upper[6] - 0.1, 0.1), 
                          m = 1, n = nStarts),
      "sigma" = comboSample(seq(bounds$lower[7] + 0.1, bounds$upper[7] - 0.1, 0.1), 
                            m = 1, n = nStarts)
    )  # We also add or subtract a little to get away from starting on the bounds
  }
  
  if (initial == TRUE){
    set.seed(101)
    parmat <- data.frame(
      "p0" = comboSample(seq(bounds$lower[1] + 1, bounds$upper[1] - 1, 2.5), 
                         m = 1, n = nStarts),
      "kg" = comboSample(seq(bounds$lower[2] + 0.1, bounds$upper[2] - 0.1, 0.1), 
                         m = 1, n = nStarts),
      "Tg" = comboSample(seq(bounds$lower[3] + 1, bounds$upper[3] - 1, 1.5), 
                         m = 1, n = nStarts),
      "kh" = comboSample(seq(bounds$lower[4] + 0.1, bounds$upper[4] - 0.1, 0.1), 
                         m = 1, n = nStarts),
      "Th" = comboSample(seq(bounds$lower[5] + 1, bounds$upper[5] - 1, 1.5), 
                         m = 1, n = nStarts),
      "Th2" = comboSample(seq(bounds$lower[6] + 0.1, bounds$upper[6] - 0.1, 0.1), 
                          m = 1, n = nStarts),
      "sigma" = comboSample(seq(bounds$lower[7] + 0.1, bounds$upper[7] - 0.1, 0.1), 
                            m = 1, n = nStarts),
      "qg" = comboSample(seq(bounds$lower[8] + 0.1, bounds$upper[8] - 0.1, 0.1), 
                         m = 1, n = nStarts),
      "qh"= comboSample(seq(bounds$lower[9] + 0.1, bounds$upper[9] - 0.1, 0.1), 
                        m = 1, n = nStarts)
    )  # We also add or subtract a little to get away from starting on the bounds
  }
  return(as.matrix(parmat))
}

# ------------------------------------------------------------------------------
# 6.33 | Function to train/test the VDR model for a given split
# ------------------------------------------------------------------------------

train_test <- function(dat, parmat, bounds, main, splits = NA, 
                       currentSplit = NA, initial){
  
  if(main == FALSE){
    # If we are training-testing on splits vs. training on the whole of block 1
    training_data <- dat[splits$train[[currentSplit]], ]
    testing_data <- dat[splits$test[[currentSplit]], ]}
  
  if(main == TRUE){
    # If training on the whole block 1 and testing on block 2
    training_data <- dat[dat$block == 1, ]
    testing_data <- dat[dat$block == 2, ]
  }
  
  # Isolate a vector of days on which measurements exist for train and test data
  measure_idx_train <- subset(training_data$day, !is.na(training_data$performance))
  measure_idx_test <- subset(testing_data$day, !is.na(testing_data$performance))
  
  # Isolate a vector of measurements in for train and test data
  measurements_train <- subset(training_data$performance, !is.na(training_data$performance))
  measurements_test <- subset(testing_data$performance, !is.na(testing_data$performance))
  
  # Put data in required format for objective function vdrObjective()
  dat_temp <- data.frame("day" = measure_idx_train, 
                         "performance" = measurements_train)
  
  load_temp <- training_data[, c("day", "load")]
  
  # Fitting iterations
  fittedModel <- optimx::multistart(parmat, 
                                    fn = vdrObjectiveLL, 
                                    lower = bounds$lower, 
                                    upper = bounds$upper,
                                    method = "L-BFGS-B",
                                    control = list(maxit = 500,
                                                   trace = FALSE),
                                    loads = load_temp,
                                    perfVals = dat_temp,
                                    initial = initial)
  
  # Compute predicted performance for the entire time-series (blocks 1 + 2)
  temp_predictions <- sapply(1:dim(parmat)[1], 
                             function(i) simulateVDR2(pars = as.numeric(fittedModel[i, 1:6]), 
                                                     loads = dat[,c("day", "load")],
                                                     if (initial == TRUE){
                                                       initialPars = as.numeric(fittedModel[i, 8:9])
                                                     } else {
                                                       initialPars = c(0,0)
                                                     }
                                                     )$performance)
  
  # Extract predictions at required days to evaluate model performance
  predictions_training <- temp_predictions[measure_idx_train, ]
  predictions_testing <- temp_predictions[measure_idx_test, ]
  
  # Compute metrics (MAPE)
  mape_train <- sapply(1:dim(predictions_training)[2],
                       function(i) mape(measured = measurements_train, 
                                        predicted = predictions_training[,i]))
  mape_test <- sapply(1:dim(predictions_testing)[2],
                      function(i) mape(measured = measurements_test,
                                       predicted = predictions_testing[,i]))
  
  # Collect these metrics and calculate average split statistics
  fittedModel <- cbind(fittedModel, mape_train, mape_test)
  
  stats <- c("mean_mape_train" = mean(mape_train), 
             "mean_mape_test" = mean(mape_test),
             "sd_mape_train" =  sd(mape_train),
             "sd_mape_test" = sd(mape_test))
  
  # Develop output
  output <- list("fittedModel" = fittedModel, "predictions" = temp_predictions,
                 "stats" = stats)
  
  return(output)
}

# ------------------------------------------------------------------------------
# 6.34 | Main CV function
# ------------------------------------------------------------------------------

expandingWindow_CV <- function(dat,
                               bounds,
                               initialWindow = 60,
                               testHorizon = 20,
                               expandRate = 4,
                               nStarts = 10,
                               cores = NULL,
                               initial = FALSE){
  
  # Generate splits (note split vectors gives you an index position vs. a 'day', think t-1!)
  splits <- generate_splits(initialWindow, testHorizon, expandRate, dat[dat$block == 1, ])
  nSplits <- length(splits$train)   # Number of splits
  
  # Create an array of random starting values over the bounds
  parmat <- create_grid(bounds, nStarts, initial = initial)
  
  # Iterate over the splits (train-test)
  if (is.null(cores)){
    cores <- detectCores(logical = TRUE)
  } # If cores not specified by user
  cl <- makeCluster(cores, type = "SOCK")              # Make cluster
  registerDoSNOW(cl)                                   # Register cluster
  fitted_splits <- foreach(i = 1:nSplits, .verbose = TRUE, .packages = c("optimx"),
                           .export = c("train_test", "vdrObjectiveLL", "simulateVDR2",
                                       "mape")) %dopar%{
                                         train_test(dat, parmat, bounds, 
                                                    main = FALSE, splits = splits,
                                                    initial = initial,
                                                    currentSplit = i)}
  # By default results returned as a list
  stopCluster(cl)                       # Stop cluster
  names(fitted_splits) <- paste0("split_",1:nSplits) # Add names to the list for each split
  
  # Compute model performance across splits and add to the existing list object
  mape_train_across <- matrix(NA, nrow = nStarts, ncol = nSplits)
  mape_test_across <- matrix(NA, nrow = nStarts, ncol = nSplits)
  for (i in 1:nSplits){
    mape_train_across[,i] <- fitted_splits[[i]]$fittedModel$mape_train
    mape_test_across[,i] <- fitted_splits[[i]]$fittedModel$mape_test
  }
  fitted_splits$across_splits <- list("training_mape" = c("mean" = mean(mape_train_across),
                                                          "sd" = sd(mape_train_across)),
                                      "testing_mape" = c("mean" = mean(mape_test_across),
                                                         "sd" = sd(mape_test_across)))
  
  fitted_splits$starting_pars <- parmat
  
  # Fit the model to all block 1 data, and test on block 2 (Main train/test split)
  mainModel <- train_test(dat, parmat, bounds, main = TRUE, initial = initial)
  
  # Extract the best set by lowest -logLik value and by prediction on test set
  mainModel$bestSet$fit <- mainModel$fittedModel[mainModel$fittedModel$value == 
                                                   min(mainModel$fittedModel$value), ]
  mainModel$bestSet$test <- mainModel$fittedModel[mainModel$fittedModel$testing_mape == 
                                                    min(mainModel$fittedModel$mape_test), ]
  
  fitted_splits$main <- mainModel
  fitted_splits$splits <- splits
  fitted_splits$nSplits <- nSplits
  fitted_splits$nStarts <- nStarts
  
  # Output results
  return(fitted_splits)
}

# ------------------------------------------------------------------------------
# 6.35 | Demonstrating the CV function call
# ------------------------------------------------------------------------------

# Call the function
example <- expandingWindow_CV(dat = dat, 
                              bounds = bounds, 
                              initialWindow = 60, 
                              testHorizon = 20,
                              expandRate = 4, 
                              nStarts = 10, 
                              cores = 1)

# Plotting and visualization analysis

  # Plotting the parameter distributions across the splits
  
    # Collate the parameters for train-test splits
    p0_across <- sapply(1:example$nSplits, function(i) rbind(example[[i]]$fittedModel$p0))
    kg_across <- sapply(1:example$nSplits, function(i) rbind(example[[i]]$fittedModel$kg))
    Tg_across <- sapply(1:example$nSplits, function(i) rbind(example[[i]]$fittedModel$Tg))
    kh_across <- sapply(1:example$nSplits, function(i) rbind(example[[i]]$fittedModel$kh))
    Th_across <- sapply(1:example$nSplits, function(i) rbind(example[[i]]$fittedModel$Th))
    Th2_across <- sapply(1:example$nSplits, function(i) rbind(example[[i]]$fittedModel$Th2))
    sigma_across <- sapply(1:example$nSplits, function(i) rbind(example[[i]]$fittedModel$sigma))
    MAPE_train_across <-  sapply(1:example$nSplits, function(i) rbind(example[[i]]$fittedModel$mape_train))
    MAPE_test_across <-  sapply(1:example$nSplits, function(i) rbind(example[[i]]$fittedModel$mape_test))
    objective_across <- sapply(1:example$nSplits, function(i) rbind(example[[i]]$fittedModel$value))
    
    # Summary of iterations on all train-test splits
    par(mfrow = c(2,5))
    boxplot(p0_across,  main = "", ylab = expression(p^{"*"}), xlab = "Split No.")
    boxplot(kg_across,  main = "", ylab = expression(k[g]), xlab = "Split No.")
    boxplot(Tg_across, main = "", ylab = expression(tau[g]), xlab = "Split No.")
    boxplot(kh_across,  main = "", ylab = expression(k[h]), xlab = "Split No.")
    boxplot(Th_across,  main = "", ylab = expression(tau[h]), xlab = "Split No.")
    boxplot(Th2_across,  main = "", ylab = expression(tau[h][2]), xlab = "Split No.")
    boxplot(sigma_across, freq = TRUE, ylab = expression(sigma), xlab = "Split No.")
    boxplot(MAPE_train_across, ylab = "MAPE (train)", xlab = "Split No.")
    boxplot(MAPE_test_across, ylab = "MAPE (test)", xlab = "Split No.")
    boxplot(objective_across, ylab = expression(-log(L)), xlab = "Split No.", main = "")
    
    # Summary of iterations on train_block1, test_block2
    par(mfrow = c(2,5))
    boxplot(example$main$fittedModel$p0, ylab = expression(p^{"*"}), main = "")
    boxplot(example$main$fittedModel$kg, ylab = expression(k[g]), main = "")
    boxplot(example$main$fittedModel$Tg, ylab = expression(tau[g]), main = "")
    boxplot(example$main$fittedModel$kh, ylab = expression(k[h]), main = "")
    boxplot(example$main$fittedModel$Th, ylab = expression(tau[h]), main = "")
    boxplot(example$main$fittedModel$Th2, ylab = expression(tau[h][2]), main = "")
    boxplot(example$main$fittedModel$sigma,  ylab = expression(sigma), main = "")
    boxplot(example$main$fittedModel$mape_train, ylab = expression("MAPE"["TRAIN"]), main = "")
    boxplot(example$main$fittedModel$mape_test, ylab = expression("MAPE"["TEST"]), main = "")
    boxplot(example$main$fittedModel$value, ylab = "-log(L)", main = "")
    
    dev.off()
    
    # Plot the main set
    plot(example$main$predictions[,1], type = "l", xlab = "Day", cex.main = 0.75,
         ylab = "Performance [a.u]", main = "Main train (block 1) - test (block 2)")
    for (i in 2:example$nStarts){
      lines(example$main$predictions[,i])
    }
    points(x = dat[dat$block == 1, "day"], y = dat[dat$block == 1, "performance"], col = "blue", pch = 1)
    points(x = dat[dat$block == 2, "day"], y = dat[dat$block == 2, "performance"], col = "red", pch = 1)
    abline(v = 100, lty = 2)
    text(x = 50, y = 95, labels = "Block 1", col = "blue")
    text(x = 128, y = 95, labels = "Block 2", col = "red")
    legend("topleft", c("Observed Data (seen)", "Observed Data (unseen)", "Model Predictions (fitted)"),
           pch = c(1, 1, NA), lty = c(NA, NA, 1), col = c("blue", "red", "black"), cex = 0.75)