library("MASS") #For mvrnorm().
m <- 1000
#Target distribution: (x, y) ~ N(mu, sigma).
mu <- c(0, 0)
rho <- 0.6
sigma <- matrix(c(1, rho, rho, 1), nrow = 2)

#Gibbs sampler: X|Y ~ N(rho*y, 1 - rho^2) and Y|X ~ N(rho*x, 1 - rho^2).

Gibbs_alg <- function() {
  X <- matrix(0, nrow = m, ncol = 2)
  X[1, ] <- c(10, 10) #Starting point.
  
  sigma <- sqrt(1 - rho^2)
  generate_x <- function(y) {
    rnorm(1, rho*y, sigma)
  }
  generate_y <- function(x) {
    rnorm(1, rho*x, sigma)
  }
  
  for(i in 1:(m - 1)) {
    X[(i + 1), 1] <- generate_x(X[i, 2])
    X[(i + 1), 2] <- generate_y(X[(i + 1), 1])
  }
  return(X)
}

X <- Gibbs_alg()

par(mfrow = c(2, 2))
plot(X[1:20, 1], X[1:20, 2], type = "l", lwd = 0.6, 
     ylim = c(-2, 10), xlim = c(-2, 10), main = "First 20 iterations")
points(X[1:20, 1], X[1:20, 2], pch = 20, lwd = 0.6, cex = 0.4)
plot(X[1:100, 1], X[1:100, 2], type = "l", lwd = 0.6, 
     ylim = c(-2, 10), xlim = c(-2, 10), main = "First 100 iterations")
points(X[1:100, 1], X[1:100, 2], pch = 20, lwd = 0.6, cex = 0.4)
plot(X[100:1000, 1], X[100:1000, 2], type = "p", pch = 20, lwd = 0.6, cex = 0.4,
     ylim = c(-3, 3), xlim = c(-3, 3), main = "Iterations 100:1000")
acf(X[100:1000, 1], 30, "correlation", T)
