# Language: R (r-project.org.uk)
# Author: Ben Stephens Hemingway
# License: GNU GPL v3
# Thesis chapter: 2
# Description: Plots in literature review (chapter 2)
# R version 4.04 (2021-02-15 "Lost Library Book")
# Package dependencies: none

# ------------------------------------------------------------------------------
# Figure 2.1 - Dose-response dynamics of the standard model for a single dose
# ------------------------------------------------------------------------------

k_g <- 1
k_h <- 1
T_g <- 30
T_h <- 10

# For right hand plot
k_g_2 <- 0.5

load <- 1

days <- c(0, 0:100)

# Banister forcing
banister_forcing_fit <- c(0, exp(-((0:100)/T_g)) * load)
banister_forcing_fat <- c(0, exp(-((0:100)/T_h)) * load)

par(mfrow = c(1,2))
# Figure 2.1-A
plot(x = days, y = banister_forcing_fit, type = "l", col = "blue",
     xlab = "Days (after dose)",
     ylab = "Dose-response level [a.u]",
     main = expression("Scaling equal;" ~ k[g] ~"="~k[h]~"="~1),
     cex.main = 0.85)
lines(x = days, y = banister_forcing_fat, col = "red", lty = 1)
lines(x = rep(T_g, 2), y = c(-1, 1/exp(1)), col = "black", lty = 2)
lines(x = rep(T_h, 2), y = c(-1, 1/exp(1)), col = "black", lty = 2)
abline(h = 1/exp(1), col = "black", lty = 2)
text(x=T_h-4, y=0, labels = expression(tau[h]))
text(x=T_g-4, y=0, labels = expression(tau[g]))
text(x=80, y=exp(-1)+0.08, labels = expression(e^-1))
legend("topright", c("Fitness", "Fatigue"), col = c("blue", "red"),
       lty = c(1,1))

# Figure 2.1-B
plot(x = days, y = banister_forcing_fat, type = "l", col = "red",
     xlab = "Days (after dose)",
     ylab = "Dose-response level [a.u]",
     main = expression("Scaling applied;" ~ k[g] ~"= 0.5,"~k[h]~"="~1),
     cex.main = 0.85)
lines(x = days, y = (banister_forcing_fit * k_g_2), col = "blue", lty = 1)
legend("topright", c("Fitness", "Fatigue"), col = c("blue", "red"),
       lty = c(1,1))

# ------------------------------------------------------------------------------
# Figure 2.3: Dose-response dynamics of the fitness-delay model for single dose
# ------------------------------------------------------------------------------

k_g <- 1
k_h <- 1
T_g1 <- 30
T_g2 <- 8
T_h <- 12

load <- 1
days <- c(0, 0:100)

calvert_forcing_fit <- c(0, (exp(-((0:100)/T_g1))-exp(-((0:100)/T_g2))) * load)
calvert_forcing_fat <- c(0,exp(-((0:100)/T_h)) * load)

gmax = max(calvert_forcing_fit)
tmax = (log(T_g1/T_g2))*((T_g1*T_g2)/(T_g1-T_g2))

par(mfrow = c(1,2))
# Figure 2.3-A
plot(x = days, y = calvert_forcing_fat, type = "l", col = "red",
     ylab = "Dose-response level (a.u)", xlab = "Days (after dose)",
     main = expression("(A)"~"Scaling equal;" ~ k[g] ~"="~k[h]~"="~1),
     cex.main = 0.85)
lines(x = days, calvert_forcing_fit, type = "l", col = "blue")
points(x = tmax, y = gmax, pch = 20, col = "black")
lines(x = rep(tmax, 2), y = c(-1, gmax), lty = 3, lwd = 0.85)
text(x=tmax + 5, y = gmax + 0.07, labels = expression("max"~g[2](t)), cex = 0.75)
text(x=tmax+1, y = 0.05, labels = expression(t[g2]~"max"), cex = 0.75)
legend("topright", c("Fitness", "Fatigue"), col = c("blue", "red"), lty = c(1,1))

# Figure 2.3-B
load <- 1/gmax # Scaling the load value for both components
calvert_forcing_fit <- c(0, (exp(-((0:100)/T_g1))-exp(-((0:100)/T_g2))) * load)
calvert_forcing_fat <- c(0,exp(-((0:100)/T_h)) * load)
calvert_forcing_fat <- calvert_forcing_fat * gmax # Scaling the fitness component
gmax2 <- max(calvert_forcing_fit)

plot(x = days, y = calvert_forcing_fat, type = "l", col = "red",
     ylab = "Dose-response level (a.u)", xlab = "Days (after dose)",
     main = expression("(B)"~omega ~"="~1~"/max"~g[2](t)~","~k[h]~"= max"~g[2](t)),
     cex.main = 0.85)
lines(x = days, calvert_forcing_fit, type = "l", col = "blue")
points(x = tmax, y = gmax2, pch = 20, col = "black")
lines(x = rep(tmax, 2), y = c(-1, gmax2), lty = 3, lwd = 0.85)
legend("topright", c("Fitness", "Fatigue"), col = c("blue", "red"), lty = c(1,1))

# ------------------------------------------------------------------------------
# Figure 2.4: Influence curves
# ------------------------------------------------------------------------------

# Parameters and day of interest (same as Fitz-Clarke et al. (1991))
T_g = 45
T_h = 15
k_g = 1
k_h = 2
t_p = 60 # Day of interest

# Function to compute the value of L(mu) over the time-series (t)

L = function(t){
  mu = t_p - t
  return((k_g * exp(-mu / T_g)) - (k_h * exp(-mu / T_h)))
}

# Calculate critical time t_n and period of maximisation t_g

t_r = t_p - (((T_g * T_h)/(T_g - T_h)) * log(k_h/k_g))
t_g = t_p - (((T_g * T_h)/(T_g - T_h)) * 
               log((k_h*T_g)/(k_g*T_h)))

# Plot the output of the function for time series 1:t_p

output <- L(1:t_p)

plot(output,type = "l", xlab = "Day", 
     ylab = expression("L" * "(" * mu * ")"), ylim = c(-1,0.5))
lines(x = rep(t_g, 2), y = c(-1, max(output)), col = "black", lty = 2)
lines(x = c(-1, t_g), y = rep(max(output),2), col = "black", lty = 2)
points(x=t_g, y = max(output), pch = 20, col = "black")
lines(x = rep(t_r, 2), y = c(-1,0), col = 'black', lty = 2)
lines(x = c(-1, t_r), y = rep(0,2), col = "black", lty = 2)
text(x=t_g - 2, y= -0.9, labels = expression(t[g]))
text(x=t_r - 2, y = -0.9, labels = expression(t[r]))
points(x=t_r, y = 0, pch = 20, col = "black")

# ------------------------------------------------------------------------------
# Figure 2.8: Hill model plots (threshold saturation)
# ------------------------------------------------------------------------------

# Changing Gamma
w_t <- seq(0,30, by = 0.2)
k <- 10
delta <- 1
gamma1 <- 5
gamma2 <- 1
gamma3 <- 0.3

hill_gamma_1 <- k * ((w_t^(gamma1))/(delta^(gamma1)+w_t^(gamma1)))
hill_gamma_2 <- k * ((w_t^(gamma2))/(delta^(gamma2)+w_t^(gamma2)))
hill_gamma_3 <- k * ((w_t^(gamma3))/(delta^(gamma3)+w_t^(gamma3)))

plot(w_t, hill_gamma_1, type = "l", main = "Hill Function", 
     xlab = "Training Load - w(t)", ylab = "Hill(w(t))")     
lines(w_t, hill_gamma_2, col = "red")    
lines(w_t, hill_gamma_3, col = "blue")

# Changing Delta
w_t <- seq(0,30, by = 0.2)
k <- 10
gamma <- 1
delta1 <- 5
delta2 <- 1
delta3 <- 0.3

hill_delta_1 <- k * ((w_t^(gamma))/(delta1^(gamma)+w_t^(gamma)))
hill_delta_2 <- k * ((w_t^(gamma))/(delta2^(gamma)+w_t^(gamma)))
hill_delta_3 <- k * ((w_t^(gamma))/(delta3^(gamma)+w_t^(gamma)))

title1 = expression(paste("Hill Function with changing ", gamma, ", ",
                          delta, " = 1, ",kappa, " = 10"))
title2 = expression(paste("Hill Function with changing ", delta, ", ",
                          gamma, " = 1, ",kappa, " = 10"))
ylab1 = expression(paste("Hill(",omega,")"))
xlab1 = expression(paste("Training Load (",omega,")"))

par(mfrow=c(1,2))
plot(w_t, hill_gamma_1, type = "l", ylim = c(0,12), col = "blue",
     main = title1, ylab = ylab1, xlab = xlab1, cex.main = 0.85, lwd = 2)     
lines(w_t, hill_gamma_2, col = "green", lty = 2, lwd = 2)    
lines(w_t, hill_gamma_3, col = "red", lty = 3, lwd = 2)
legend("bottomright",c(
  expression(paste(gamma," = 5")),
  expression(paste(gamma," = 1")),
  expression(paste(gamma," = 0.3"))
), lty = 1:3, col = c("blue", "green", "red"), lwd = rep(2,3))
plot(w_t, hill_delta_1, type = "l", ylim = c(0,12), col = "blue", lwd = 2,
     main = title2, ylab = ylab1, xlab = xlab1, cex.main = 0.85)     
lines(w_t, hill_delta_2, col = "green", lty = 2, lwd = 2)    
lines(w_t, hill_delta_3, col = "red", lty = 3, lwd = 2)
legend("bottomright",c(
  expression(paste(delta," = 5")),
  expression(paste(delta," = 1")),
  expression(paste(delta," = 0.3"))
), lty = 1:3, col = c("blue", "green", "red"), lwd = rep(2,3))