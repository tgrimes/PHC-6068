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
  #X[1, ] <- c(0, 0) #Starting value.
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


#############################
# Random-walk 
#############################
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




#############################
# Independence Sampler
#############################
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



#############################
# Langevin Diffusion Process
#############################
#Proposal distribution T(X -> Y):
#The Langevin diffusion process: y = x + s^2/2 del(log(pi(x))) + s*epsilon.
sigma <- 1.4 #Step size.
#Gradient of log(f(x))
del_logpi <- function(X) {
  x <- X[1]
  y <- X[2]
  f = pi(X)
  del_x = (1/f)*(-4*x*exp(-0.5*(4*x^2 + 0.5*y^2)) - 
                    8*(x - 5)*exp(-0.5*(4*(x - 5)^2 + 0.5*(y - 5)^2)))
  del_y = (1/f)*(-0.5*y*exp(-0.5*(4*x^2 + 0.5*y^2)) - 
                    (y - 5)*exp(-0.5*(4*(x - 5)^2 + 0.5*(y - 5)^2)))
  return(c(del_x, del_y))
}
#Density of the proposal (up to a proportion).
q <- function(X, Y) {
  del <- del_logpi(X)
  mean <- X + sigma^2/2*del
  exp(-(t(Y - mean) %*% (Y - mean))/(2*sigma^2))
}
#Generator for the proposal distribution.
rq <- function(n = 1, X) {
  del <- del_logpi(X)
  epsilon <- mvrnorm(1, c(0, 0), diag(2))
  return(X + sigma^2/2*del + sigma*epsilon)
}

X3 <- MH_alg(q, rq, m)



#Plotting the results:
par(mfrow = c(2, 3))
plot_xy <- function(X, title) {
  plot(X[, 1], X[, 2], type = "l", lwd = 0.05, ylim = c(-4, 10), xlim = c(-2, 8),
       main = title, xlab = "x", ylab = "y")
  points(X[, 1], X[, 2], pch = 20, lwd = 0.1, cex = 0.1)
}
plot_x <- function(X) { 
  plot(1:m, X[, 1], type = "l", lwd = 0.4, ylim = c(-2, 8),
       xlab = "iteration", ylab = "x") 
}

plot_xy(X1, "Random-walk")
plot_xy(X2, "Independence Sampler")
plot_xy(X3, "Langevin diffusion")

plot_x(X1)
plot_x(X2)
plot_x(X3)