library("MASS") #For mvrnorm().
m <- 20000

#Density of interest (up to a proportion).
pi <- function(X) {
  x <- X[1]
  y <- X[2]
  exp(-0.5*(4*x^2 + 0.5*y^2)) + 2*exp(-0.5*(4*(x - 5)^2 + 0.5*(y - 5)^2))
}

#############################
# M-H algorithm.
#############################
MH_alg <- function(q, rq, m = 20000) {
  #M-H ratio:
  ratio <- function(X, Y) {
    pi(Y)*q(Y, X)/(pi(X)*q(X, Y))
  }
  
  X <- matrix(0, nrow = m, ncol = 2)
  X[1, ] <- c(0, 0) #Starting value.
  for(i in 1:(m - 1)) {
    Y <- rq(1, X[i, ])
    if(runif(1) < min(1, ratio(X[i, ], Y))) {
      X[i + 1, ] <- Y
    } else {
      X[i + 1, ] <- X[i, ]
    }
  }
  
  return(X)
}

#Proposal distribution T(X -> Y):
#We will use an random-walk proposal T(X, Y) = q(Y) with Y|X ~ N(X, sigma)
sigma <- 2*diag(2)
#Density of the proposal (up to a proportion).
q <- function(X, Y) {
  Y <- Y 
  mean <- X
  exp(-0.5*(t(Y - mean) %*% solve(sigma^2) %*% (Y - mean)))
}
#Generator for the proposal distribution.
rq <- function(n = 1, X) {
  mean <- X
  mvrnorm(n, mean, sigma)
}

X1 <- MH_alg(q, rq, m)



#Proposal distribution T(X -> Y):
#We will use an independence sampler T(X, Y) = q(Y), with Y|X ~ N(mean, sigma)
mean <- c(2.5, 2.5)
sigma <- 4*diag(2)
#Density of the proposal (up to a proportion).
q <- function(X, Y) {
  exp(-0.5*(t(Y - mean) %*% solve(sigma^2) %*% (Y - mean)))
}
#Generator for the proposal distribution.
rq <- function(n = 1, X) {
  mvrnorm(n, mean, sigma)
}

X2 <- MH_alg(q, rq, m)

par(mfrow = c(2, 2))
plot(X1[, 1], X1[, 2], type = "l", lwd = 0.1, ylim = c(-4, 10), xlim = c(-2, 8),
     main = "Random-walk", xlab = "x", ylab = "y")
points(X1[, 1], X1[, 2], pch = 20, lwd = 0.1, cex = 0.1)
plot(X2[, 1], X2[, 2], type = "l", lwd = 0.1, ylim = c(-4, 10), xlim = c(-2, 8),
     main = "Independence sampler", xlab = "x", ylab = "y")
points(X2[, 1], X2[, 2], pch = 20, lwd = 0.1, cex = 0.1)

plot(1:m, X1[, 1], type = "l", lwd = 0.4, ylim = c(-2, 8),
     xlab = "iteration", ylab = "x")
plot(1:m, X2[, 1], type = "l", lwd = 0.4, ylim = c(-2, 8),
     xlab = "iteration", ylab = "x")

