# Admin

  source("analysis_fn.R")
  library(plyr)
  dir.create("analysis_output")
  dir.create("analysis_output/plots")

# -----------------------------
# Standard Model
# -----------------------------  
  
# Load data
  
  s1 <- readRDS("raw/standard/1/gridFull.Rda")
  s2 <- readRDS("raw/standard/2/gridFull.Rda")
  s3 <- readRDS("raw/standard/3/gridFull.Rda")

  s1 <- as.data.frame(s1)
  s2 <- as.data.frame(s2)
  s3 <- as.data.frame(s3)
 
# ------------------------------------------------------------------------------
# STANDARD MODEL SCENARIOS
# ------------------------------------------------------------------------------

# Separate the Data Sets

  # Separate sets that converged to true parameters (based on differences with true pars)
  s1_true <- s1[s1$p0TrueDiff < 0.01 & s1$kgTrueDiff < 0.01 & s1$TgTrueDiff < 0.1 & 
                s1$khTrueDiff < 0.01 & s1$ThTrueDiff < 0.1 & s1$convCode == 0, ]
  s2_true <- s2[s2$p0TrueDiff < 0.01 & s2$kgTrueDiff < 0.01 & s2$TgTrueDiff < 0.1 & 
                s2$khTrueDiff < 0.01 & s2$ThTrueDiff < 0.1 & s2$convCode == 0, ]
  s3_true <- s3[s3$p0TrueDiff < 0.01 & s3$kgTrueDiff < 0.01 & s3$TgTrueDiff < 0.1 &
                s3$khTrueDiff < 0.01 & s3$ThTrueDiff < 0.1 & s3$convCode == 0, ]
  
  # Separate the sets that had abnormal convergence/termination in L-BFGS-B
  s1_error <- s1[s1$convCode != 0,]
  s2_error <- s2[s2$convCode != 0,]
  s3_error <- s3[s3$convCode != 0,]
  
  # Separate sets that converged (according to L-BFGS-B) at other points
  s1_alt <- s1[(s1$p0TrueDiff < 0.1 & s1$kgTrueDiff < 0.01 & s1$TgTrueDiff < 0.1 & 
                s1$khTrueDiff < 0.01 & s1$ThTrueDiff < 0.1) == FALSE & s1$convCode == 0, ]
  s2_alt <- s2[(s2$p0TrueDiff < 0.1 & s2$kgTrueDiff < 0.01 & s2$TgTrueDiff < 0.1 & 
                s2$khTrueDiff < 0.01 & s2$ThTrueDiff < 0.1) == FALSE & s2$convCode == 0, ]
  s3_alt <- s3[(s3$p0TrueDiff < 0.1 & s3$kgTrueDiff < 0.01 & s3$TgTrueDiff < 0.1 & 
                s3$khTrueDiff < 0.01 & s3$ThTrueDiff < 0.1) == FALSE & s3$convCode == 0, ]
   
# Distributional statistics on the obtained separated sets (min, max, median, mad of each parameter)
  
  # SCENARIO 1 (100% Data, Standard Model)
  # Sets that converged to true parameters - check min/max equal and m.a.d. < 0.001
  s1_true_stats <- cbind(apply(s1_true[,6:12], 2, min), apply(s1_true[,6:12], 2, max), 
                         apply(s1_true[,6:12], 2, median), apply(s1_true[,6:12], 2, mad))
  colnames(s1_true_stats) <- c("min","max","median","mad")
  
  # Sets that converged to other local extrema
  s1_alt_stats <- cbind(apply(s1_alt[,6:12], 2, min), apply(s1_alt[,6:12], 2, max), 
                        apply(s1_alt[,6:12], 2, median), apply(s1_alt[,6:12], 2, mad))
  colnames(s1_alt_stats) <- c("min","max","median","mad")
  
  # Sets where the L-BFGS-B algorithm terminated unsuccessfully according to its convergence criteria
  s1_error_stats <- cbind(apply(s1_error[,6:12], 2, min), apply(s1_error[,6:12], 2, max), 
                          apply(s1_error[,6:12], 2, median), apply(s1_error[,6:12], 2, mad))
  colnames(s1_error_stats) <- c("min","max","median","mad")
  
  # SCENARIO 2 (50% Data, Standard Model)
  # Sets that converged to true parameters - check min/max equal and m.a.d. < 0.001
  s2_true_stats <- cbind(apply(s2_true[,6:12], 2, min), apply(s2_true[,6:12], 2, max), 
                         apply(s2_true[,6:12], 2, median), apply(s2_true[,6:12], 2, mad))
  colnames(s2_true_stats) <- c("min","max","median","mad")
  
  # Sets that converged to other local extrema
  s2_alt_stats <- cbind(apply(s2_alt[,6:12], 2, min), apply(s2_alt[,6:12], 2, max), 
                        apply(s2_alt[,6:12], 2, median), apply(s2_alt[,6:12], 2, mad))
  colnames(s2_alt_stats) <- c("min","max","median","mad")
  
  # Sets where the L-BFGS-B algorithm terminated unsuccessfully according to its convergence criteria
  s2_error_stats <- cbind(apply(s2_error[,6:12], 2, min), apply(s2_error[,6:12], 2, max), 
                          apply(s2_error[,6:12], 2, median), apply(s2_error[,6:12], 2, mad))
  colnames(s2_error_stats) <- c("min","max","median","mad")
  
  # SCENARIO 3 (33% Data, Standard Model)
  # -----------------------------------------------------------------------------------------
  # Sets that converged to true parameters - check min/max equal and m.a.d. < 0.001
  s3_true_stats <- cbind(apply(s3_true[,6:12], 2, min), apply(s3_true[,6:12], 2, max), 
                         apply(s3_true[,6:12], 2, median), apply(s3_true[,6:12], 2, mad))
  colnames(s3_true_stats) <- c("min","max","median","mad")
  
  # Sets that converged to other local extrema
  s3_alt_stats <- cbind(apply(s3_alt[,6:12], 2, min), apply(s3_alt[,6:12], 2, max), 
                        apply(s3_alt[,6:12], 2, median), apply(s3_alt[,6:12], 2, mad))
  colnames(s3_alt_stats) <- c("min","max","median","mad")
  
  # Sets where the L-BFGS-B algorithm terminated unsuccessfully according to its convergence criteria
  s3_error_stats <- cbind(apply(s3_error[,6:12], 2, min), apply(s3_error[,6:12], 2, max), 
                          apply(s3_error[,6:12], 2, median), apply(s3_error[,6:12], 2, mad))
  colnames(s3_error_stats) <- c("min","max","median","mad")
  
  # Combine the results from scenario 1,2,3 and write to output directory
  # -----------------------------------------------------------------------------------------
  s1_stats <- rbind(t(s1_true_stats),t(s1_alt_stats),t(s1_error_stats))
  s2_stats <- rbind(t(s2_true_stats),t(s2_alt_stats),t(s2_error_stats))
  s3_stats <- rbind(t(s3_true_stats),t(s3_alt_stats),t(s3_error_stats))
  s_stats <- rbind(s1_stats,s2_stats,s3_stats)
  s_stats[,1] <- round(s_stats[,1], digits = 0)
  s_stats[,2:6] <- round(s_stats[,2:6],digits = 3)
  s_stats[,7] <- round(s_stats[,7], digits = 2)
  s_stats <- as.data.frame(s_stats)
  s_stats$scenario <- c(rep("s1", 12), rep("s2",12), rep("s3",12))
  s_stats$set <- rep(c(rep("true",4),rep("alt",4),rep("error",4)),3)
  write.csv(s_stats, "analysis_output/standard_parameter_distributions.csv")

# Evaluate Critical points (i.e., local minima or saddle)
  
  # s1 (100% fitting data)
  s1_true$type <- ifelse (s1_true$EigenVal1 > 0 & s1_true$EigenVal2 > 0 & 
                            s1_true$EigenVal3 > 0 & s1_true$EigenVal4 > 0 & 
                            s1_true$EigenVal5 > 0, "minimum", 
                          ifelse(s1_true$EigenVal1 < 0 | s1_true$EigenVal2 < 0 | 
                                 s1_true$EigenVal3 < 0 | s1_true$EigenVal4 < 0 | 
                                 s1_true$EigenVal5 < 0, "saddle", "unknown"))
  s1_alt$type <- ifelse (s1_alt$EigenVal1 > 0 & s1_alt$EigenVal2 > 0 & s1_alt$EigenVal3 > 0 & 
                         s1_alt$EigenVal4 > 0 & s1_alt$EigenVal5 > 0, "minimum", 
                            ifelse(s1_alt$EigenVal1 < 0 | s1_alt$EigenVal2 < 0 | 
                                   s1_alt$EigenVal3 < 0 | s1_alt$EigenVal4 < 0 | 
                                   s1_alt$EigenVal5 < 0, "saddle", "unknown"))
  
  # s2 (50% fitting data)
  s2_true$type <- ifelse (s2_true$EigenVal1 > 0 & s2_true$EigenVal2 > 0 & 
                          s2_true$EigenVal3 > 0 & s2_true$EigenVal4 > 0 & 
                          s2_true$EigenVal5 > 0, "minimum", 
                          ifelse(s2_true$EigenVal1 < 0 | s2_true$EigenVal2 < 0 | 
                          s2_true$EigenVal3 < 0 | s2_true$EigenVal4 < 0 | 
                          s2_true$EigenVal5 < 0, "saddle", "unknown"))
  s2_alt$type <- ifelse (s2_alt$EigenVal1 > 0 & s2_alt$EigenVal2 > 0 & 
                         s2_alt$EigenVal3 > 0 & s2_alt$EigenVal4 > 0 & 
                         s2_alt$EigenVal5 > 0, "minimum", 
                           ifelse(s2_alt$EigenVal1 < 0 | s2_alt$EigenVal2 < 0 | 
                                  s2_alt$EigenVal3 < 0 | s2_alt$EigenVal4 < 0 | 
                                  s2_alt$EigenVal5 < 0, "saddle", "unknown"))
  
  # s3 (33% fitting data)
  s3_true$type <- ifelse (s3_true$EigenVal1 > 0 & s3_true$EigenVal2 > 0 & 
                          s3_true$EigenVal3 > 0 & s3_true$EigenVal4 > 0 & 
                          s3_true$EigenVal5 > 0, "minimum", 
                            ifelse(s3_true$EigenVal1 < 0 | s3_true$EigenVal2 < 0 | 
                                   s3_true$EigenVal3 < 0 | s3_true$EigenVal4 < 0 | 
                                   s3_true$EigenVal5 < 0, "saddle", "unknown"))
  s3_alt$type <- ifelse (s3_alt$EigenVal1 > 0 & s3_alt$EigenVal2 > 0 & 
                         s3_alt$EigenVal3 > 0 & s3_alt$EigenVal4 > 0 & 
                         s3_alt$EigenVal5 > 0, "minimum", 
                           ifelse(s3_alt$EigenVal1 < 0 | s3_alt$EigenVal2 < 0 | 
                                  s3_alt$EigenVal3 < 0 | s3_alt$EigenVal4 < 0 | 
                                  s3_alt$EigenVal5 < 0, "saddle", "unknown"))
  
  # Compile the results of this sorting process and write to output directory
  critical_points <- data.frame("model" = c("standard", 'standard', "standard"),
                                "data" = c("100%", "50%", "33%"),
                                "data_points" = c(147,74,49),
                                "iterations" = c(10^5, 10^5, 10^5),
                                "true" = c(sum(s1_true$type == "minimum"), 
                                           sum(s2_true$type == "minimum"), 
                                           sum(s3_true$type == "minimum")),
                                "local_min" = c(sum(s1_alt$type == "minimum"), 
                                                sum(s2_alt$type == "minimum"), 
                                                sum(s3_alt$type == "minimum")),
                                "saddle" = c(sum(s1_alt$type == "saddle"), 
                                             sum(s2_alt$type == "saddle"), 
                                             sum(s3_alt$type == "saddle")))
  write.csv(critical_points, "analysis_output/standard_critical_points_types.csv")

# Identify uniqueness in the sets among the alternative critical points found (tables)
  
  # Standard model (100%)
  s1_unique <- s1_alt
  s1_unique$p0Fitted <- round(s1_unique$p0Fitted, digits = 1)
  s1_unique$kgFitted <- round(s1_unique$kgFitted, digits = 2)
  s1_unique$TgFitted <- round(s1_unique$TgFitted, digits = 1)
  s1_unique$khFitted <- round(s1_unique$khFitted, digits = 2)
  s1_unique$ThFitted <- round(s1_unique$ThFitted, digits = 1)
  s1_unique <- count(s1_unique, vars = c("p0Fitted","kgFitted", 
                                         "TgFitted","khFitted",
                                         "ThFitted","type"))
  s1_unique <- s1_unique[order(s1_unique$freq, decreasing = TRUE),]
  s1_unique$rss <- apply(s1_unique[,1:5], 1, standardLoss, 
                         standardCriterion = standardCriterion1, 
                         loads = loads, standardIndex = standardIndex1)
  write.csv(s1_unique, "analysis_output/standard_1_unique_sets.csv")
  
  # Standard model (50%)
  s2_unique <- s2_alt
  s2_unique$p0Fitted <- round(s2_unique$p0Fitted, digits = 1)
  s2_unique$kgFitted <- round(s2_unique$kgFitted, digits = 2)
  s2_unique$TgFitted <- round(s2_unique$TgFitted, digits = 1)
  s2_unique$khFitted <- round(s2_unique$khFitted, digits = 2)
  s2_unique$ThFitted <- round(s2_unique$ThFitted, digits = 1)
  s2_unique <- count(s2_unique, vars = c("p0Fitted","kgFitted",
                                         "TgFitted","khFitted","ThFitted","type"))
  s2_unique <- s2_unique[order(s2_unique$freq, decreasing = TRUE),]
  s2_unique$rss <- apply(s2_unique[,1:5], 1, standardLoss, 
                         standardCriterion = standardCriterion2, 
                         loads = loads, standardIndex = standardIndex2)
  write.csv(s2_unique, "analysis_output/standard_2_unique_sets.csv")
  
  # Standard model (33%)
  s3_unique <- s2_alt
  s3_unique$p0Fitted <- round(s3_unique$p0Fitted, digits = 1)
  s3_unique$kgFitted <- round(s3_unique$kgFitted, digits = 2)
  s3_unique$TgFitted <- round(s3_unique$TgFitted, digits = 1)
  s3_unique$khFitted <- round(s3_unique$khFitted, digits = 2)
  s3_unique$ThFitted <- round(s3_unique$ThFitted, digits = 1)
  s3_unique <- count(s3_unique, vars = c("p0Fitted","kgFitted",
                                         "TgFitted","khFitted",
                                         "ThFitted","type"))
  s3_unique <- s3_unique[order(s3_unique$freq, decreasing = TRUE),]
  s3_unique$rss <- apply(s3_unique[,1:5], 1, standardLoss, 
                         standardCriterion = standardCriterion3, 
                         loads = loads, standardIndex = standardIndex3)
  write.csv(s3_unique, "analysis_output/standard_3_unique_sets.csv")
  

# Tabulate the prediction errors (write to output directory)
  
  # Scenario 1 (33% data)
  s1_predictions <- rbind(
    t(cbind(apply(s1_true[,22:23], 2, min), apply(s1_true[,22:23], 2, max), 
            apply(s1_true[,22:23], 2, median), apply(s1_true[,22:23], 2, mad))),
    t(cbind(apply(s1_alt[,22:23], 2, min), apply(s1_alt[,22:23], 2, max), 
            apply(s1_alt[,22:23], 2, median), apply(s1_alt[,22:23], 2, mad))),
    t(cbind(apply(s1_error[,22:23], 2, min), apply(s1_error[,22:23], 2, max), 
            apply(s1_error[,22:23], 2, median), apply(s1_error[,22:23], 2, mad)))
  )
  colnames(s1_predictions) <- c("RMSE","MAPE")
  rownames(s1_predictions) <- rep(c("min","max","median","mad"),3)
  
  # Scenario 2 (50% data)
  s2_predictions <- rbind(
    t(cbind(apply(s2_true[,22:23], 2, min), apply(s2_true[,22:23], 2, max), 
            apply(s2_true[,22:23], 2, median), apply(s2_true[,22:23], 2, mad))),
    t(cbind(apply(s2_alt[,22:23], 2, min), apply(s2_alt[,22:23], 2, max), 
            apply(s2_alt[,22:23], 2, median), apply(s2_alt[,22:23], 2, mad))),
    t(cbind(apply(s2_error[,22:23], 2, min), apply(s2_error[,22:23], 2, max), 
            apply(s2_error[,22:23], 2, median), apply(s2_error[,22:23], 2, mad)))
  )
  colnames(s2_predictions) <- c("RMSE","MAPE")
  rownames(s2_predictions) <- rep(c("min","max","median","mad"),3)
  
  s3_predictions <- rbind(
    t(cbind(apply(s3_true[,22:23], 2, min), apply(s3_true[,22:23], 2, max), 
            apply(s3_true[,22:23], 2, median), apply(s3_true[,22:23], 2, mad))),
    t(cbind(apply(s3_alt[,22:23], 2, min), apply(s3_alt[,22:23], 2, max), 
            apply(s3_alt[,22:23], 2, median), apply(s3_alt[,22:23], 2, mad))),
    t(cbind(apply(s3_error[,22:23], 2, min), apply(s3_error[,22:23], 2, max), 
            apply(s3_error[,22:23], 2, median), apply(s3_error[,22:23], 2, mad)))
  )
  colnames(s3_predictions) <- c("RMSE","MAPE")
  rownames(s3_predictions) <- rep(c("min","max","median","mad"),3)
  standard_predictions <- as.data.frame(rbind(s1_predictions, s2_predictions, s3_predictions))
  standard_predictions$data <- c(rep("100",12),rep("50",12),rep("33",12))
  standard_predictions$set <- rep(c(rep("true",4),rep("alternate",4),rep("error",4)),3)
  standard_predictions[,1:2] <- round(standard_predictions[,1:2], digits = 3)
  write.csv(standard_predictions, "analysis_output/standard_model_fit.csv")
  
# Plotting the predictions
  
  # Boxplots (RMSE, MAPE) - Standard Model all three scenarios
  par(mfrow = c(1,2))
  boxplot(list("100% Data" = s1_alt[,22], "50% Data" = s2_alt[,22], 
               "33% Data" = s3_alt[,22]), main = "RMSE")
  boxplot(list("100% Data" = s1_alt[,23], "50% Data" = s2_alt[,23], 
               "33% Data" = s3_alt[,23]), main = "MAPE")
  
  # Plot the RSS (cost function) distributions
  par(mfrow = c(1,2))
  boxplot(list("100% Data" = s1_alt$optimRSS, "50% Data" = s2_alt$optimRSS, 
               "33% Data" = s3_alt$optimRSS), 
          main = expression("RSS - Standard Model"), 
          ylab = "RSS (cost)", cex.main = 0.85)
  boxplot(list("100% Data" = s1_alt[s1_alt$optimRSS < 500, "optimRSS"], 
               "50% Data" = s2_alt[s2_alt$optimRSS < 500, "optimRSS"], 
               "33% Data" = s3_alt[s3_alt$optimRSS < 500, "optimRSS"]), 
          main = "RSS - Standard Model (zoomed in, < 500)", 
          ylab = "RSS (cost)", cex.main = 0.85)
  
  
  # Compute the predictions (alternative solutions only)
  s1_performance <- apply(s1_alt[,7:11], 1, standardCompute, 
                          loads = loads, returnObject = "performance")
  s2_performance <- apply(s2_alt[,7:11], 1, standardCompute, 
                          loads = loads, returnObject = "performance")
  s3_performance <- apply(s3_alt[,7:11], 1, standardCompute, 
                          loads = loads, returnObject = "performance")
  
  # Plot prediction traces
  
  par(mfrow  = c(3,1))
  matplot(s1_performance, xlab = "Day", ylab = "Performance (a.u)",
          main = "Scenario: 100% Fitting Data x Standard Model | Fitted (non-true) solutions",
          cex.main = 0.9, type = "l", col = c("chartreuse1"))
  points(standardCriterion1, pch = 1, cex = 0.8)
  legend("bottomright", c("Fitting data", "Fitted model predictions"), pch = c(1,NA), 
         lty = c(NA,1), lwd = c(NA,2), col = c("black","chartreuse1"),
         cex = 0.95, bty = "n")
  matplot(s2_performance, xlab = "Day", ylab = "Performance (a.u)",
          main = "Scenario: 50% Fitting Data x Standard Model | Fitted (non-true) solutions",
          cex.main = 0.9, type = "l", col = c("chartreuse1"))
  points(standardIndex2, standardCriterion2, pch = 1, cex = 0.8)
  legend("bottomright", c("Fitting data", "Fitted model predictions"), pch = c(1,NA), 
         lty = c(NA,1), lwd = c(NA,2), col = c("black","chartreuse1"),
         cex = 0.95, bty = "n")
  matplot(s3_performance, xlab = "Day", ylab = "Performance (a.u)",
          main = "Scenario: 33% Fitting Data x Standard Model | Fitted (non-true) solutions",
          cex.main = 0.9, type = "l", col = c("chartreuse1"))
  points(standardIndex3, standardCriterion3, pch = 1, cex = 0.8)
  legend("bottomright", c("Fitting data", "Fitted model predictions"), pch = c(1,NA), 
         lty = c(NA,1), lwd = c(NA,2), col = c("black","chartreuse1"),
         cex = 0.95, bty = "n")


# ------------------------------------------------------------------------------
# CALVERT MODEL SCENARIOS (FITNESS-DELAY MODEL)
# ------------------------------------------------------------------------------
  
# Load Data
  
  c1 <- readRDS("raw/calvert/1/gridFull.Rda")
  c2 <- readRDS("raw/calvert/2/gridFull.Rda")
  c3 <- readRDS("raw/calvert/3/gridFull.Rda")
  c1 <- as.data.frame(c1)
  c2 <- as.data.frame(c2)
  c3 <- as.data.frame(c3)
  
# Segment the data-sets
  
  # Separate sets that converged to true parameters (based on differences with true pars)
  c1_true <- c1[c1$p0TrueDiff < 0.1 & c1$kgTrueDiff < 0.1 & c1$Tg1TrueDiff < 0.1 & 
                c1$Tg2TrueDiff < 0.1 & c1$khTrueDiff < 0.1 & c1$ThTrueDiff < 0.1 & 
                c1$convCode == 0, ]
  c2_true <- c2[c2$p0TrueDiff < 0.1 & c2$kgTrueDiff < 0.1 & c2$Tg1TrueDiff < 0.1 & 
                c2$Tg2TrueDiff < 0.1 & c2$khTrueDiff < 0.1 & c2$ThTrueDiff < 0.1 & 
                c2$convCode == 0, ]
  c3_true <- c3[c3$p0TrueDiff < 0.1 & c3$kgTrueDiff < 0.1 & c3$Tg1TrueDiff < 0.1 & 
                c3$Tg2TrueDiff < 0.1 & c3$khTrueDiff < 0.1 & c3$ThTrueDiff < 0.1 & 
                c3$convCode == 0, ]
  
  # Separate sets that had abnormal convergence/termination
  c1_error <- c1[c1$convCode != 0, ]
  c2_error <- c2[c2$convCode != 0, ]
  c3_error <- c3[c3$convCode != 0, ]
  
  # Separate sets that converged or terminated at other points
  c1_alt <- c1[(c1$p0TrueDiff < 0.1 & c1$kgTrueDiff < 0.1 & c1$Tg1TrueDiff < 0.1 & 
                c1$Tg2TrueDiff < 0.1 & c1$khTrueDiff < 0.1 & c1$ThTrueDiff < 0.1) == FALSE & c1$convCode == 0, ]
  c2_alt <- c1[(c2$p0TrueDiff < 0.1 & c2$kgTrueDiff < 0.1 & c2$Tg1TrueDiff < 0.1 & 
                c2$Tg2TrueDiff < 0.1 & c2$khTrueDiff < 0.1 & c2$ThTrueDiff < 0.1) == FALSE & c2$convCode == 0, ]
  c3_alt <- c3[(c3$p0TrueDiff < 0.1 & c3$kgTrueDiff < 0.1 & c3$Tg1TrueDiff < 0.1 & 
                c3$Tg2TrueDiff < 0.1 & c3$khTrueDiff < 0.1 & c3$ThTrueDiff < 0.1) == FALSE & c3$convCode == 0, ]
  
# Distributional statistics on the obtained separated sets (min, max, median, mad of each parameter)
  
  c1_true_stats <- cbind(apply(c1_true[,7:14], 2, min), apply(c1_true[,7:14], 2, max), 
                         apply(c1_true[,7:14], 2, median), apply(c1_true[,7:14], 2, mad))
  colnames(c1_true_stats) <- c("min","max","median","mad")
  c1_alt_stats <- cbind(apply(c1_alt[,7:14], 2, min), apply(c1_alt[,7:14], 2, max), 
                        apply(c1_alt[,7:14], 2, median), apply(c1_alt[,7:14], 2, mad))
  colnames(c1_alt_stats) <- c("min","max","median","mad")
  c1_error_stats <- cbind(apply(c1_error[,7:14], 2, min), apply(c1_error[,7:14], 2, max), 
                          apply(c1_error[,7:14], 2, median), apply(c1_error[,7:14], 2, mad))
  colnames(c1_error_stats) <- c("min","max","median","mad")
  c2_true_stats <- cbind(apply(c2_true[,7:14], 2, min), apply(c2_true[,7:14], 2, max), 
                         apply(c2_true[,7:14], 2, median), apply(c2_true[,7:14], 2, mad))
  colnames(c2_true_stats) <- c("min","max","median","mad")
  c2_alt_stats <- cbind(apply(c2_alt[,7:14], 2, min), apply(c2_alt[,7:14], 2, max), 
                        apply(c2_alt[,7:14], 2, median), apply(c2_alt[,7:14], 2, mad))
  colnames(c2_alt_stats) <- c("min","max","median","mad")
  c2_error_stats <- cbind(apply(c2_error[,7:14], 2, min), apply(c2_error[,7:14], 2, max), 
                          apply(c2_error[,7:14], 2, median), apply(c2_error[,7:14], 2, mad))
  colnames(c2_error_stats) <- c("min","max","median","mad")
  c3_true_stats <- cbind(apply(c3_true[,7:14], 2, min), apply(c3_true[,7:14], 2, max), 
                         apply(c3_true[,7:14], 2, median), apply(c3_true[,7:14], 2, mad))
  colnames(c3_true_stats) <- c("min","max","median","mad")
  c3_alt_stats <- cbind(apply(c3_alt[,7:14], 2, min), apply(c3_alt[,7:14], 2, max), 
                        apply(c3_alt[,7:14], 2, median), apply(c3_alt[,7:14], 2, mad))
  colnames(c3_alt_stats) <- c("min","max","median","mad")
  c3_error_stats <- cbind(apply(c3_error[,7:14], 2, min), apply(c3_error[,7:14], 2, max), 
                          apply(c3_error[,7:14], 2, median), apply(c3_error[,7:14], 2, mad))
  colnames(c3_error_stats) <- c("min","max","median","mad")
  
  # Output
  c1_stats <- rbind(t(c1_true_stats),t(c1_alt_stats),t(c1_error_stats))
  c2_stats <- rbind(t(c2_true_stats),t(c2_alt_stats),t(c2_error_stats))
  c3_stats <- rbind(t(c3_true_stats),t(c3_alt_stats),t(c3_error_stats))
  c_stats <- rbind(c1_stats,c2_stats,c3_stats)
  c_stats[,1] <- round(c_stats[,1], digits = 0)
  c_stats[,2:6] <- round(c_stats[,2:6],digits = 3)
  c_stats[,7] <- round(c_stats[,7], digits = 2)
  c_stats <- as.data.frame(c_stats)
  c_stats$scenario <- c(rep("s1", 12), rep("s2",12), rep("s3",12))
  c_stats$set <- rep(c(rep("true",4),rep("alt",4),rep("error",4)),3)
  write.csv(c_stats, "analysis_output/calvert_parameter_distributions.csv")
  
# Evaluate Critical points (local minima or saddle)
  
  # Calvert (100%)
  c1_true$type <- ifelse (c1_true$EigenVal1 > 0 & c1_true$EigenVal2 > 0 & 
                            c1_true$EigenVal3 > 0 & c1_true$EigenVal4 > 0 & 
                            c1_true$EigenVal5 > 0 & c1_true$EigenVal6 > 0, "minimum", 
                              ifelse(c1_true$EigenVal1 < 0 | c1_true$EigenVal2 < 0 | 
                                       c1_true$EigenVal3 < 0 | c1_true$EigenVal4 < 0 | 
                                       c1_true$EigenVal5 < 0 | c1_true$EigenVal6 < 0, 
                                     "saddle", "unknown"))
  c1_alt$type <- ifelse (c1_alt$EigenVal1 > 0 & c1_alt$EigenVal2 > 0 & 
                           c1_alt$EigenVal3 > 0 & c1_alt$EigenVal4 > 0 & 
                           c1_alt$EigenVal5 > 0 & c1_alt$EigenVal6 > 0, "minimum", 
                              ifelse(c1_alt$EigenVal1 < 0 | c1_alt$EigenVal2 < 0 | 
                                     c1_alt$EigenVal3 < 0 | c1_alt$EigenVal4 < 0 | 
                                     c1_alt$EigenVal5 < 0 | c1_alt$EigenVal6 < 0, 
                                     "saddle", "unknown"))
  
  # Calvert (50%)
  c2_true$type <- ifelse (c2_true$EigenVal1 > 0 & c2_true$EigenVal2 > 0 & 
                          c2_true$EigenVal3 > 0 & c2_true$EigenVal4 > 0 & 
                          c2_true$EigenVal5 > 0 & c2_true$EigenVal6 > 0, "minimum", 
                            ifelse(c2_true$EigenVal1 < 0 | c2_true$EigenVal2 < 0 | 
                                    c2_true$EigenVal3 < 0 | c2_true$EigenVal4 < 0 | 
                                    c2_true$EigenVal5 < 0 | c2_true$EigenVal6 < 0, 
                                   "saddle", "unknown"))
  c2_alt$type <- ifelse (c2_alt$EigenVal1 > 0 & c2_alt$EigenVal2 > 0 & 
                          c2_alt$EigenVal3 > 0 & c2_alt$EigenVal4 > 0 & 
                          c2_alt$EigenVal5 > 0 & c2_alt$EigenVal6 > 0, "minimum", 
                         ifelse(c2_alt$EigenVal1 < 0 | c2_alt$EigenVal2 < 0 | 
                                c2_alt$EigenVal3 < 0 | c2_alt$EigenVal4 < 0 | 
                                  c2_alt$EigenVal5 < 0 | c2_alt$EigenVal6 < 0, 
                                "saddle", "unknown"))
  
  # Calvert (33)
  c3_true$type <- ifelse (c3_true$EigenVal1 > 0 & c3_true$EigenVal2 > 0 & 
                          c3_true$EigenVal3 > 0 & c3_true$EigenVal4 > 0 & 
                          c3_true$EigenVal5 > 0 & c3_true$EigenVal6 > 0, "minimum", 
                              ifelse(c3_true$EigenVal1 < 0 | c3_true$EigenVal2 < 0 | 
                                     c3_true$EigenVal3 < 0 | c3_true$EigenVal4 < 0 | 
                                     c3_true$EigenVal5 < 0 | c3_true$EigenVal6 < 0, 
                                     "saddle", "unknown"))
  c3_alt$type <- ifelse (c3_alt$EigenVal1 > 0 & c3_alt$EigenVal2 > 0 & 
                         c3_alt$EigenVal3 > 0 & c3_alt$EigenVal4 > 0 & 
                         c3_alt$EigenVal5 > 0 & c3_alt$EigenVal6 > 0, "minimum", 
                         ifelse(c3_alt$EigenVal1 < 0 | c3_alt$EigenVal2 < 0 | 
                                c3_alt$EigenVal3 < 0 | c3_alt$EigenVal4 < 0 | 
                                c3_alt$EigenVal5 < 0 | c3_alt$EigenVal6 < 0, 
                                "saddle", "unknown"))
  
  # Output
  critical_points <- data.frame("model" = c("calvert", 'calvert', "calvert"),
                                "data" = c("100%", "50%", "33%"),
                                "data_points" = c(147,74,49),
                                "iterations" = c(7^6, 7^6, 7^6),
                                "true" = c(sum(c1_true$type == "minimum"), 
                                           sum(c2_true$type == "minimum"), 
                                           sum(c3_true$type == "minimum")),
                                "local_min" = c(sum(c1_alt$type == "minimum"), 
                                                sum(c2_alt$type == "minimum"), 
                                                sum(c3_alt$type == "minimum")),
                                "saddle" = c(sum(c1_alt$type == "saddle"), 
                                             sum(c2_alt$type == "saddle"), 
                                             sum(c3_alt$type == "saddle")))
  write.csv(critical_points, "analysis_output/calvert_critical_points_types.csv")

# Identify unique sets among the alternative critical points (tables)
  
  # Standard model (100%)
  c1_unique <- c1_alt
  c1_unique$p0Fitted <- round(c1_unique$p0Fitted, digits = 1)
  c1_unique$kgFitted <- round(c1_unique$kgFitted, digits = 2)
  c1_unique$Tg1Fitted <- round(c1_unique$Tg1Fitted, digits = 1)
  c1_unique$Tg2Fitted <- round(c1_unique$Tg2Fitted, digits = 1)
  c1_unique$khFitted <- round(c1_unique$khFitted, digits = 2)
  c1_unique$ThFitted <- round(c1_unique$ThFitted, digits = 1)
  c1_unique <- count(c1_unique, vars = c("p0Fitted","kgFitted",
                                         "Tg1Fitted","Tg2Fitted",
                                         "khFitted","ThFitted","type"))
  c1_unique <- c1_unique[order(c1_unique$freq, decreasing = TRUE),]
  c1_unique$rss <- apply(c1_unique[,1:6], 1, calvertLoss, 
                         calvertCriterion = calvertCriterion1, 
                         loads = loads, calvertIndex = calvertIndex1)
  write.csv(c1_unique, "analysis_output/calvert_1_unique_sets.csv")
  
  # Standard model (50%)
  c2_unique <- c2_alt
  c2_unique$p0Fitted <- round(c2_unique$p0Fitted, digits = 1)
  c2_unique$kgFitted <- round(c2_unique$kgFitted, digits = 2)
  c2_unique$Tg1Fitted <- round(c2_unique$Tg1Fitted, digits = 1)
  c2_unique$Tg2Fitted <- round(c2_unique$Tg2Fitted, digits = 1)
  c2_unique$khFitted <- round(c2_unique$khFitted, digits = 2)
  c2_unique$ThFitted <- round(c2_unique$ThFitted, digits = 1)
  c2_unique <- count(c2_unique, vars = c("p0Fitted","kgFitted",
                                         "Tg1Fitted","Tg2Fitted",
                                         "khFitted","ThFitted","type"))
  c2_unique <- c2_unique[order(c2_unique$freq, decreasing = TRUE),]
  c2_unique$rss <- apply(c2_unique[,1:6], 1, calvertLoss, 
                         calvertCriterion = calvertCriterion2, 
                         loads = loads, calvertIndex = calvertIndex2)
  write.csv(c2_unique, "analysis_output/calvert_2_unique_sets.csv")
  
  # Standard model (33%)
  c3_unique <- c3_alt
  c3_unique$p0Fitted <- round(c3_unique$p0Fitted, digits = 1)
  c3_unique$kgFitted <- round(c3_unique$kgFitted, digits = 2)
  c3_unique$Tg1Fitted <- round(c3_unique$Tg1Fitted, digits = 1)
  c3_unique$Tg2Fitted <- round(c3_unique$Tg2Fitted, digits = 1)
  c3_unique$khFitted <- round(c3_unique$khFitted, digits = 2)
  c3_unique$ThFitted <- round(c3_unique$ThFitted, digits = 1)
  c3_unique <- count(c3_unique, vars = c("p0Fitted","kgFitted",
                                         "Tg1Fitted","Tg2Fitted",
                                         "khFitted","ThFitted","type"))
  c3_unique <- c3_unique[order(c3_unique$freq, decreasing = TRUE),]
  c3_unique$rss <- apply(c3_unique[,1:6], 1, calvertLoss, 
                         calvertCriterion = calvertCriterion3, 
                         loads = loads, calvertIndex = calvertIndex3)
  write.csv(c3_unique, "analysis_output/calvert_3_unique_sets.csv")

# Plot the RSS distributions
  
  par(mfrow = c(1,2))
  boxplot(list("100% Data" = c1_alt$optimRSS, "50% Data" = c2_alt$optimRSS, 
               "33% Data" = c3_alt$optimRSS), 
          main = expression("RSS - Fitness-delay Model"), 
          ylab = "RSS (cost)", cex.main = 0.85)
  boxplot(list("100% Data" = c1_alt[c1_alt$optimRSS < 500, "optimRSS"], 
               "50% Data" = c2_alt[c2_alt$optimRSS < 500, "optimRSS"], 
               "33% Data" = c3_alt[c3_alt$optimRSS < 500, "optimRSS"]), 
          main = "RSS - Fitness-delay Model (zoomed in, < 500)", 
          ylab = "RSS (cost)", cex.main = 0.85)
  
# Tabulate the prediction errors
  
  c1_predictions <- rbind(
    t(cbind(apply(c1_true[,25:26], 2, min), apply(c1_true[,25:26], 2, max), 
            apply(c1_true[,25:26], 2, median), apply(c1_true[,25:26], 2, mad))),
    t(cbind(apply(c1_alt[,25:26], 2, min), apply(c1_alt[,25:26], 2, max), 
            apply(c1_alt[,25:26], 2, median), apply(c1_alt[,25:26], 2, mad))),
    t(cbind(apply(c1_error[,25:26], 2, min), apply(c1_error[,25:26], 2, max), 
            apply(c1_error[,25:26], 2, median), apply(c1_error[,25:26], 2, mad)))
  )
  colnames(c1_predictions) <- c("RMSE","MAPE")
  rownames(c1_predictions) <- rep(c("min","max","median","mad"),3)
  
  c2_predictions <- rbind(
    t(cbind(apply(c2_true[,25:26], 2, min), apply(c2_true[,25:26], 2, max), 
            apply(c2_true[,25:26], 2, median), apply(c2_true[,25:26], 2, mad))),
    t(cbind(apply(c2_alt[,25:26], 2, min), apply(c2_alt[,25:26], 2, max), 
            apply(c2_alt[,25:26], 2, median), apply(c2_alt[,25:26], 2, mad))),
    t(cbind(apply(c2_error[,25:26], 2, min), apply(c2_error[,25:26], 2, max), 
            apply(c2_error[,25:26], 2, median), apply(c2_error[,25:26], 2, mad)))
  )
  colnames(c2_predictions) <- c("RMSE","MAPE")
  rownames(c2_predictions) <- rep(c("min","max","median","mad"),3)
  
  c3_predictions <- rbind(
    t(cbind(apply(c3_true[,25:26], 2, min), apply(c3_true[,25:26], 2, max),
            apply(c3_true[,25:26], 2, median), apply(c3_true[,25:26], 2, mad))),
    t(cbind(apply(c3_alt[,25:26], 2, min), apply(c3_alt[,25:26], 2, max), 
            apply(c3_alt[,25:26], 2, median), apply(c3_alt[,25:26], 2, mad))),
    t(cbind(apply(c3_error[,25:26], 2, min), apply(c3_error[,25:26], 2, max), 
            apply(c3_error[,25:26], 2, median), apply(c3_error[,25:26], 2, mad)))
  )
  colnames(c3_predictions) <- c("RMSE","MAPE")
  rownames(c3_predictions) <- rep(c("min","max","median","mad"),3)
  
  calvert_predictions <- as.data.frame(rbind(c1_predictions, c2_predictions, c3_predictions))
  calvert_predictions$data <- c(rep("100",12),rep("50",12),rep("33",12))
  calvert_predictions$set <- rep(c(rep("true",4),rep("alternate",4),rep("error",4)),3)
  calvert_predictions[,1:2] <- round(calvert_predictions[,1:2], digits = 3)
  write.csv(calvert_predictions, "analysis_output/calvert_model_fit.csv")
  
  # Model fit boxplots
  par(mfrow = c(1,2))
  boxplot(list("100% Data" = c1_alt[,25], "50% Data" = c2_alt[,25], 
               "33% Data" = c3_alt[,25]), main = "RMSE")
  boxplot(list("100% Data" = c1_alt[,26], "50% Data" = c2_alt[,26], 
               "33% Data" = c3_alt[,26]), main = "MAPE")
  
  # Compute the predictions (alternative solutions only)
  c1_performance <- apply(c1_alt[,8:13], 1, calvertCompute, 
                          loads = loads, returnObject = "performance")
  c2_performance <- apply(c2_alt[,8:13], 1, calvertCompute, 
                          loads = loads, returnObject = "performance")
  c3_performance <- apply(c3_alt[,8:13], 1, calvertCompute, 
                          loads = loads, returnObject = "performance")

  # Plotting the prediction traces
  par(mfrow  = c(3,1))
  matplot(c1_performance, xlab = "Day", ylab = "Performance (a.u)",
          main = "Scenario: 100% Fitting Data x Fitness-delay Model | Fitted (non-true) solutions",
          cex.main = 0.95, type = "l", col = c("orange"))
  points(calvertCriterion1, pch = 1, cex = 0.9)
  legend("bottomright", c("Fitting data", "Fitted model predictions"), pch = c(1,NA), 
         lty = c(NA,1), lwd = c(NA,2), col = c("black","orange"),
         cex = 0.95, bty = "n")
  matplot(c2_performance, xlab = "Day", ylab = "Performance (a.u)",
          main = "Scenario: 50% Fitting Data x Fitness-delay Model | Fitted (non-true) solutions",
          cex.main = 0.95, type = "l", col = c("orange"))
  points(calvertIndex2, calvertCriterion2, pch = 1, cex = 0.9)
  legend("bottomright", c("Fitting data", "Fitted model predictions"), pch = c(1,NA), 
         lty = c(NA,1), lwd = c(NA,2), col = c("black","orange"),
         cex = 0.95, bty = "n")
  matplot(c3_performance, xlab = "Day", ylab = "Performance (a.u)",
          main = "Scenario: 33% Fitting Data x Fitness-delay Model | Fitted (non-true) solutions",
          cex.main = 0.95, type = "l", col = c("orange"))
  points(calvertIndex3, calvertCriterion3, pch = 1, cex = 0.9)
  legend("bottomright", c("Fitting data", "Fitted model predictions"), pch = c(1,NA), 
         lty = c(NA,1), lwd = c(NA,2), col = c("black","orange"),
         cex = 0.95, bty = "n")
  
# ------------------------------------------------------------------------------
# COMBINED PLOT - BOTH MODEL SCENARIOS (PARAMETERS)  
# ------------------------------------------------------------------------------  

  # Plot the parameter distributions
  
  par(mfrow = c(6,2), mar = c(2.5,2.5,2.5,2.5))
  boxplot(list("100% Data" = s1_alt$p0Fitted, "50% Data" = s2_alt$p0Fitted, 
               "33% Data" = s3_alt$p0Fitted), 
          main = expression(p~'*' ~' - Standard Model'), cex.main = 1)
  abline(h=100, col = "red")
  boxplot(list("100% Data" = s1_alt$kgFitted, "50% Data" = s2_alt$kgFitted, 
               "33% Data" = s3_alt$kgFitted), 
          main = expression(k[g]~' - Standard Model'), cex.main = 1)
  abline(h=0.72, col = "red")
  boxplot(list("100% Data" = s1_alt$TgFitted, "50% Data" = s2_alt$TgFitted, 
               "33% Data" = s3_alt$TgFitted), 
          main = expression(tau[g]~' - Standard Model'), cex.main = 1)
  abline(h=28.5, col = "red")
  boxplot(list("100% Data" = s1_alt$khFitted, "50% Data" = s2_alt$khFitted, 
               "33% Data" = s3_alt$khFitted), 
          main = expression(k[h]~' - Standard Model'), cex.main = 1)
  abline(h=1.2, col = "red")
  boxplot(list("100% Data" = s1_alt$ThFitted, "50% Data" = s2_alt$ThFitted, 
               "33% Data" = s3_alt$ThFitted), 
          main = expression(tau[h]~' - Standard Model'), cex.main = 1)
  abline(h=8.6, col = "red")
  boxplot(list("100% Data" = c1_alt$p0Fitted, "50% Data" = c2_alt$p0Fitted, 
               "33% Data" = c3_alt$p0Fitted), 
          main = expression(p~'*' ~' - Fitness-delay Model'), cex.main = 1)
  abline(h=100, col = "red")
  boxplot(list("100% Data" = c1_alt$kgFitted, "50% Data" = c2_alt$kgFitted, 
               "33% Data" = c3_alt$kgFitted), 
          main = expression(k[g]~' - Fitness-delay Model'), cex.main = 1)
  abline(h=0.72, col = "red")
  boxplot(list("100% Data" = c1_alt$Tg1Fitted, "50% Data" = c2_alt$Tg1Fitted, 
               "33% Data" = c3_alt$Tg1Fitted), 
          main = expression(tau[g][1]~' - Fitness-delay Model'), cex.main = 1)
  abline(h=32.5, col = "red")
  boxplot(list("100% Data" = c1_alt$Tg2Fitted, "50% Data" = c2_alt$Tg2Fitted, 
               "33% Data" = c3_alt$Tg2Fitted), 
          main = expression(tau[g][2]~' - Fitness-delay Model'), cex.main = 1)
  abline(h=4.3, col = "red")
  boxplot(list("100% Data" = c1_alt$khFitted, "50% Data" = c2_alt$khFitted, 
               "33% Data" = c3_alt$khFitted), 
          main = expression(k[h]~' - Fitness-delay Model'), cex.main = 1)
  abline(h=1.05, col = "red")
  boxplot(list("100% Data" = c1_alt$ThFitted, "50% Data" = c2_alt$ThFitted, 
               "33% Data" = c3_alt$ThFitted), 
          main = expression(tau[h]~' - Fitness-delay Model'), cex.main = 1)
  abline(h=8.6, col = "red")
  
  
# Initial distributions (Appendices)
  
  par(mfrow = c(3,5))
  boxplot(s1_true[,1], main = "p*")
  boxplot(s1_true[,2], main = expression(k[g]))
  boxplot(s1_true[,3], main = expression(tau[g]))
  boxplot(s1_true[,4], main = expression(k[h]))
  boxplot(s1_true[,5], main = expression(tau[h]))
  boxplot(s2_true[,1], main = "p*")
  boxplot(s2_true[,2], main = expression(k[g]))
  boxplot(s2_true[,3], main = expression(tau[g]))
  boxplot(s2_true[,4], main = expression(k[h]))
  boxplot(s2_true[,5], main = expression(tau[h]))
  boxplot(s3_true[,1], main = "p*")
  boxplot(s3_true[,2], main = expression(k[g]))
  boxplot(s3_true[,3], main = expression(tau[g]))
  boxplot(s3_true[,4], main = expression(k[h]))
  boxplot(s3_true[,5], main = expression(tau[h]))
  
  par(mfrow = c(1,3))
  boxplot(s1_true[,6], main = "Init RSS (100% Data)")
  boxplot(s2_true[,6], main = "Init RSS (50% Data)")
  boxplot(s3_true[,6], main = "Init RSS (33% Data")
  
  par(mfrow = c(3,6))
  boxplot(c1_true[,1], main = "p*")
  boxplot(c1_true[,2], main = expression(k[g]))
  boxplot(c1_true[,3], main = expression(tau[g[1]]))
  boxplot(c1_true[,4], main = expression(tau[g[1]]))
  boxplot(c1_true[,5], main = expression(k[h]))
  boxplot(c1_true[,6], main = expression(tau[h]))
  boxplot(c2_true[,1], main = "p*")
  boxplot(c2_true[,2], main = expression(k[g]))
  boxplot(c2_true[,3], main = expression(tau[g[1]]))
  boxplot(c2_true[,4], main = expression(tau[g[1]]))
  boxplot(c2_true[,5], main = expression(k[h]))
  boxplot(c2_true[,6], main = expression(tau[h]))
  boxplot(c3_true[,1], main = "p*")
  boxplot(c3_true[,2], main = expression(k[g]))
  boxplot(c3_true[,3], main = expression(tau[g[1]]))
  boxplot(c3_true[,4], main = expression(tau[g[1]]))
  boxplot(c3_true[,5], main = expression(k[h]))
  boxplot(c3_true[,6], main = expression(tau[h]))
  
  par(mfrow = c(1,3))
  boxplot(c1_true[,6], main = "Init RSS (100% Data)")
  boxplot(c2_true[,6], main = "Init RSS (50% Data)")
  boxplot(c3_true[,6], main = "Init RSS (33% Data")
  
  
  