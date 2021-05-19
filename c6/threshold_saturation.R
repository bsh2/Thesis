# ******************************************************************************
# Language: R (r-project.org.uk)
# Author: Ben Stephens Hemingway
# License: GNU GPL v3
# Thesis chapter: 6
# Subsection: 6.3
# Description: Threshold saturation - plots
# R version 4.04 (2021-02-15 "Lost Library Book")
# Package dependencies: None
# ******************************************************************************

# Introduce the Hill function
hillTransform <- function(kappa, delta, gamma, loads){
  hill_w <- kappa * (loads^(gamma) / (delta^(gamma) + loads^(gamma)))
  return(hill_w)
}

# Develop plots to illustrate behavior
plot(x=seq(0,3,0.1),
     hillTransform(3, 1, 3, loads = seq(0,3,0.1)), type = "l", col = "red",
     xlab = expression(omega), ylab = expression("Hill("~omega~")"), lwd = 1.5,
     main = expression(kappa ~"="~3~"(fixed)"))
lines(x=seq(0,3,0.1),
      hillTransform(3,1,2, loads = seq(0,3,0.1)), type = "l", col = "green",
      lty = 2, lwd = 1.5)
lines(x=seq(0,3,0.1),
      hillTransform(3,1,1, loads = seq(0,3,0.1)), type = "l", col = "blue",
      lty = 3, lwd = 1.5)
legend("topleft", c(expression(delta ~"="~1~"," ~ gamma ~"="~3),
                    expression(delta ~"="~1~"," ~ gamma ~"="~2),
                    expression(delta ~"="~1~"," ~ gamma ~"="~1)),
       lty = c(1,2,3), lwd = 1.5, col = c("red", "green", "blue"), cex = 0.85)

plot(x=seq(0,3,0.1),
     hillTransform(3, 3, 1, loads = seq(0,3,0.1)), type = "l", col = "red",
     xlab = expression(omega), ylab = expression("Hill("~omega~")"), lwd = 1.5,
     main = expression(kappa ~"="~3~"(fixed)"), ylim = c(0,3.5))
lines(x=seq(0,3,0.1),
      hillTransform(3,2,1, loads = seq(0,3,0.1)), type = "l", col = "green",
      lty = 2, lwd = 1.5)
lines(x=seq(0,3,0.1),
      hillTransform(3,0.5,1, loads = seq(0,3,0.1)), type = "l", col = "blue",
      lty = 3, lwd = 1.5)
legend("topleft", c(expression(delta ~"="~3~"," ~ gamma ~"="~1),
                    expression(delta ~"="~2~"," ~ gamma ~"="~1),
                    expression(delta ~"="~0.5~"," ~ gamma ~"="~1)),
       lty = c(1,2,3), lwd = 1.5, col = c("red", "green", "blue"), cex = 0.85)

par(mfrow = c(1,2))
plot(x=seq(0,3,0.1),
     hillTransform(3, 1, 6, loads = seq(0,3,0.1)), type = "l", col = "red",
     xlab = expression(omega), ylab = expression("Hill("~omega~")"), lwd = 1.5,
     main = expression(kappa ~"="~3~"(fixed)"), ylim = c(0,3.5), lty = 1)
lines(x=seq(0,3,0.1),
      hillTransform(3,1,5, loads = seq(0,3,0.1)), type = "l", col = "green",
      lty = 2, lwd = 1.5)
lines(x=seq(0,3,0.1),
      hillTransform(3,1,4, loads = seq(0,3,0.1)), type = "l", col = "orange",
      lty = 2, lwd = 1.5)
lines(x=seq(0,3,0.1),
      hillTransform(3,1,3, loads = seq(0,3,0.1)), type = "l", col = "blue",
      lty = 2, lwd = 1.5)
lines(x=seq(0,3,0.1),
      hillTransform(3,1,2, loads = seq(0,3,0.1)), type = "l", col = "pink",
      lty = 2, lwd = 1.5)
lines(x=seq(0,3,0.1),
      hillTransform(3,1,1, loads = seq(0,3,0.1)), type = "l", col = "purple",
      lty = 2, lwd = 1.5)
legend("topleft", c(expression(delta ~"="~1~"," ~ gamma ~"="~6),
                    expression(delta ~"="~1~"," ~ gamma ~"="~5),
                    expression(delta ~"="~1~"," ~ gamma ~"="~4),
                    expression(delta ~"="~1~"," ~ gamma ~"="~3),
                    expression(delta ~"="~1~"," ~ gamma ~"="~2),
                    expression(delta ~"="~1~"," ~ gamma ~"="~1)),
       lty = rep(1,6), lwd = 1.5, col = c("red", "green", "orange",
                                          "blue", "pink", "purple"), cex = 0.7)
plot(x=seq(0,3,0.1),
     hillTransform(3, 2, 6, loads = seq(0,3,0.1)), type = "l", col = "red",
     xlab = expression(omega), ylab = expression("Hill("~omega~")"), lwd = 1.5,
     main = expression(kappa ~"="~3~"(fixed)"), ylim = c(0,3.5), lty = 1)
lines(x=seq(0,3,0.1),
      hillTransform(3,2,5, loads = seq(0,3,0.1)), type = "l", col = "green",
      lty = 2, lwd = 1.5)
lines(x=seq(0,3,0.1),
      hillTransform(3,2,4, loads = seq(0,3,0.1)), type = "l", col = "orange",
      lty = 2, lwd = 1.5)
lines(x=seq(0,3,0.1),
      hillTransform(3,2,3, loads = seq(0,3,0.1)), type = "l", col = "blue",
      lty = 2, lwd = 1.5)
lines(x=seq(0,3,0.1),
      hillTransform(3,2,2, loads = seq(0,3,0.1)), type = "l", col = "pink",
      lty = 2, lwd = 1.5)
lines(x=seq(0,3,0.1),
      hillTransform(3,2,1, loads = seq(0,3,0.1)), type = "l", col = "purple",
      lty = 2, lwd = 1.5)
legend("topleft", c(expression(delta ~"="~2~"," ~ gamma ~"="~6),
                    expression(delta ~"="~2~"," ~ gamma ~"="~5),
                    expression(delta ~"="~2~"," ~ gamma ~"="~4),
                    expression(delta ~"="~2~"," ~ gamma ~"="~3),
                    expression(delta ~"="~2~"," ~ gamma ~"="~2),
                    expression(delta ~"="~2~"," ~ gamma ~"="~1)),
       lty = rep(1,6), lwd = 1.5, col = c("red", "green", "orange",
                                          "blue", "pink", "purple"), cex = 0.75)