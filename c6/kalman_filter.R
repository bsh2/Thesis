# ******************************************************************************
# Language: R (r-project.org.uk)
# Author: Ben Stephens Hemingway & BA-Ogorek
# License: GNU GPL v3
# Thesis chapter: 6
# Subsection: 6.4
# Description: Kalman filtering (state-space reformulation)
# R version 4.04 (2021-02-15 "Lost Library Book")
# Package dependencies:
library(MASS)
library(GA)
# ******************************************************************************

# ------------------------------------------------------------------------------
# 6.37 | Defining a function for instantiating a state-space FFM
# ------------------------------------------------------------------------------

state_space_FFM <- function(pars){
  
  # Transition matrix (eq. 6.37)
  A <- matrix(c(exp(-1 / pars$tau_g), 0, 0, exp(-1 / pars$tau_h)), ncol = 2)
  
  # State intercept (eq. 6.38)
  B <- matrix(c(exp(-1 / pars$tau_g), exp(-1 / pars$tau_h)), ncol = 1)
  
  # Measurement matrix
  C <- matrix(c(pars$k_g, -1 * pars$k_h), ncol = 2)
  
  # Variances (eq's 6.36, 6.35)
  Q <- matrix(c(pars$sigma_g^2, rep(pars$rho_gh * pars$sigma_g * pars$sigma_h, 2), 
                pars$sigma_h^2), ncol = 2)
  xi <- pars$xi
  
  # Prior distribution of fitness and fatigue (initial conditions)
  x_0 <- c(pars$g_0, pars$h_0)
  M_0 <- matrix(c(pars$sd_g0^2, 
                  rep(pars$rho_gh0 * pars$sd_g0 * pars$sd_h0, 2), pars$sd_h0^2), 
                ncol = 2)
  
  model <- list(A = A, B = B, C = C, Q = Q, xi = xi, x_0 = x_0,
                M_0 = M_0, p_star = pars$p_star)
  
  return(model)
}

# ------------------------------------------------------------------------------
# 6.38 | Defining a function to simulate the state-space model
# ------------------------------------------------------------------------------

simulate_ss_model <- function(ss_model, loads){
  
  # Set up vectors and state matrix structures
  T <- length(loads$load)
  performance <- numeric(length(loads$load))
  X <- matrix(rep(NA, 2 * T), ncol = 2)
  
  # Simulate (note %*% is matrix multiplication operator in R)
  for (n in 1:T){
    # A priori mean and variance of state - (eq.6.38)
    if (n == 1){
      X[n, ] <- ss_model$x_0  # Unconditional: x_0
    } else{
      # Conditional: x_n | x_(n-1) - (eq. 6.34)
      X[n, ] <- mvrnorm(1, (ss_model$A %*% X[n - 1, ] + 
                            ss_model$B * loads$load[n - 1]), ss_model$Q)
      # Uses MASS package function mvrnorm(n, mu, sigma) to draw from
    }
    # Simulate conditional system state: p_n | x_n - (eq. 6.39)
    performance[n] <- rnorm(1, ss_model$p_star + ss_model$C %*% X[n, ], ss_model$xi)
  }
  
  simulation <- data.frame("day" = 1:T, "load" = loads$load, "performance" = performance,
                           "fitness" = X[,1], "fatigue" = X[,2])
  return(simulation)
}

# ------------------------------------------------------------------------------
# 6.39 | Simulating the state-space FFM under mock data with random noise terms
# ------------------------------------------------------------------------------

loads <- rep(c(1, 1.2, 0.5, 1.8, 2, 0.25, 0.7, 0.9, 0, 0.5, 1, 0.8, 1.2, 1.3, 
               3, 0.9, 0, 0, 2, 1.1), 5)
loads <- data.frame("day" = 0:length(loads), "loads" = c(0,loads))

pars <- list(p_star = 95, k_g = 0.85, k_h = 1.2, tau_g = 26, tau_h = 5, g_0 = 10,
             h_0 = 6, xi = 2.5, sd_g0 = 2, sd_h0 = 1, rho_gh0 = 0, sigma_g = 2,
             sigma_h = 1.2, rho_gh = 0.3)

# Instantiate Kalman model and simulate under the parameters
ss_model <- state_space_FFM(pars)
set.seed(109)
simulated_ss_model <- simulate_ss_model(ss_model, loads)
plot(simulated_ss_model$performance, type = "l", ylab = "Performance [a.u]",
     xlab = "Day", main = "Simulated state-space FFM with random noise")

# ------------------------------------------------------------------------------
# 6.40 | Kalman filtering function
# ------------------------------------------------------------------------------

kalman_filter <- function(ss_model, dat){
  
  # Extract properties of the dataset
  T <- nrow(dat)
  loads <- dat$load
  measured <- dat$performance
  
  # Set up object structures
  loglike <- numeric(T)
  X <- matrix(rep(NA, 2 * T), ncol = 2)  # a posteriori state estimate
  M <- matrix(rep(NA, 4 * T), ncol = 4)  # Vectorised state vcovs
  Z <- matrix(rep(NA, 2 * T), ncol = 2)  # a priori state estimates
  
  # Kalman updating equations
  for (n in 1:T){
    if (n == 1){   # Initialisation
      z_n <- ss_model$A %*% ss_model$x_0 + ss_model$B * mean(loads)
      P_n <- ss_model$Q + ss_model$A %*% ss_model$M_0 %*% t(ss_model$A)
    } else{
      z_n <- ss_model$A %*% X[n-1, ] + ss_model$B * loads[n-1]
      P_n <- ss_model$Q + ss_model$A %*% matrix(M[n-1,], ncol = 2) %*% t(ss_model$A) 
    }
    
    # Likelihood of performance measurement
    S_n <- ss_model$xi^2 + ss_model$C %*% P_n %*% t(ss_model$C)
    e_n <- measured[n] - (ss_model$p_star + ss_model$C %*% z_n)
    loglike[n] <- dnorm(e_n, mean = 0, sd = sqrt(S_n), log = TRUE)
    
    # A posterori mean and variance
    K_n <- P_n %*% t(ss_model$C) %*% (1 / S_n) # Kalman gain
    X[n, ] <- z_n + K_n %*% e_n
    M[n, ] <- as.vector((diag(2) - K_n %*% ss_model$C) %*% P_n)
    Z[n, ] <- z_n
    
  } # End for loop
  p_hat <- ss_model$p_star + X %*% t(ss_model$C) # Filtered predictions of performance
  output <- list(p_hat = p_hat, g_hat = X[, 1], h_hat = X[, 2], M = M, Z = Z,
                 loglike = loglike)
  return(output)
}

# ------------------------------------------------------------------------------
# 6.41 | Demonstration of Kalman filtering for the simulated data
# ------------------------------------------------------------------------------

filtered_model <- kalman_filter(ss_model, dat = simulated_ss_model)

# Plot the results
plot(filtered_model$p_hat, type = "l", lty = 3, col = "blue", ylab = "Performance [a.u]",
     xlab = "Day", lwd = 2)
points(simulated_ss_model$performance, lty = 1, pch = 1, col = "red")
legend("topleft", c("Kalman filtered FFM", "Simulated data"),
       lty = c(3, NA), pch = c(NA, 1), col = c("blue", "red"), lwd = c(3, NA), cex = 0.75)

# ------------------------------------------------------------------------------
# 6.42 | Defining a function to fit a state-space model with Kalman filter to data
# ------------------------------------------------------------------------------

fit_filtered_model <- function(dat, box){
  
  extract_pars <- function(par){
    pars <- list(p_star = par[1],
                 k_g = par[2], 
                 k_h = par[3],
                 tau_g = par[4],
                 tau_h = par[5],
                 g_0 = par[6],
                 h_0 = par[7],
                 xi = par[8],
                 sd_g0 = par[9],
                 sd_h0 = par[10],
                 rho_gh0 = 0, # Fix rho_gh0 to zero
                 sigma_g = par[11],
                 sigma_h = par[12],
                 rho_gh = par[13])
    return(pars)
  }
  
  log_likelihood <- function(par){
    pars <- extract_pars(par)
    temp_model <- state_space_FFM(pars)
    filtered_model <- kalman_filter(ss_model = temp_model, dat = dat)
    # Note as GA is a maximiser by default, we just need the log likelihood
    return(sum(filtered_model$loglike))
  }
  
  maximise_likelihood <- GA::ga("real-valued",
                                fitness = log_likelihood,
                                lower = box$lower,
                                upper = box$upper,
                                maxiter = 2000,
                                monitor = TRUE,
                                optim = TRUE, # Local search via L-BFGS-B (stochastic)
                                optimArgs = list(method = "L-BFGS-B",
                                                 lower = box$lower,
                                                 upper = box$upper,
                                                 poptim = 0.2, # Probability of LS
                                                 pressel = 0.3),
                                popSize = 100,
                                parallel = TRUE
                                )

  pars <- extract_pars(maximise_likelihood@solution)
  plot(maximise_likelihood)
  
  fitted_model <- state_space_FFM(pars)
  return(list("fitted_model" = fitted_model, "optimisation" = maximise_likelihood))
}



# c(p*, kg, kh, Tau_g, Tau_h, g_0, h_0, xi, sd_g0, sd_h0, sigma_g, sigma_h, rho_gh)
box_constraints <- data.frame("lower" = c(50, 0.1, 0.1, 1, 1, 1, 0.5, 0.5, 0.01, 0.01,
                              1, 1, -0.999), "upper" = c(150, 5, 5, 50, 50, 20, 
                                                         20, 10, 5, 5, 5, 5, 0.999))

fitted_kalman_model <- fit_filtered_model(dat = simulated_ss_model, box = box_constraints)
fitted_kalman_predictions <- kalman_filter(fitted_kalman_model$fitted_model, simulated_ss_model)

# Output (fitted parameters, predicted values, plots)
print(fitted_kalman_model)
print(fitted_kalman_predictions)

plot(simulated_ss_model$day, simulated_ss_model$performance, pch = 1, col = "red",
     ylab = "Performance [a.u]", xlab = "Day")
lines(fitted_kalman_predictions$p_hat, type = "l", col = "blue", lty = 3, lwd = 1.5)
legend("topleft", c("Measured data", "Kalman model (fitted)"),
       lty = c(NA, 3), lwd = c(NA, 1.5), pch = c(1, NA), col = c("red", "blue"))


