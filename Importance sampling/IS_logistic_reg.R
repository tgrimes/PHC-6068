library("MASS") #For mvrnorm().


seed <- floor(10^4 * runif(1))
set.seed(seed) #Good: 5147  Bad: 886 (right, close), 1650 (left), 4138 (right)
n <- 50 #Sample size.
m <- 20000 #IS simulation size.
k <- 500 #Number of times to resample.
beta <- c(1, 2)
p <- length(beta)
sigma <- 100 #Standard deviation of proposed dist for beta.



#############################
# Generating a random sample
#############################
generate_sample <- function(beta) {
  #X's are iid standard normal.
  x <- cbind(rep(1, n), mvrnorm(n, rep(0, (p - 1)), diag(p - 1)))
  
  #Conditional distribution of Y|X and beta.
  prob_y <- function(x, beta) {
    prob <- 1/(exp(-(x %*% beta)) + 1)
    return(prob)
  }
  y <- apply(x, 1, function(x) { rbinom(1, 1, prob_y(x, beta)) })
  # n <- nrow(x)
  # #Generate n uniform random variables.
  # u <- runif(n, 0, 1)
  # #Calculate each conditional probability of y given x and beta.
  # prob <- prob_y(x, beta)
  # 
  # #Use the inverse-cdf method to generate y values. 
  # y <- rep(0, n)
  # y[u < prob] <- 1
  
  return(cbind(y, x))
}

####################################
# Importance Sampling procedure:
####################################
IS_estimates <- function(y, x) {
  #Sample beta:
  beta_sampled <- mvrnorm(m, rep(0, p), sigma*diag(p)) 
  
  prob_y <- 1 + exp(-(x %*% t(beta_sampled))) #n by m.
  prob_y[y == 0, ] <- (1 + exp(x %*% t(beta_sampled)))[y == 0, ]
  prob_y <- apply(prob_y, 2, prod) 
  
  #Calculate the individual weights and sum of weights.
  weight <- 1/prob_y
  sum_weight <- sum(weight)
  
  E_beta_hat <- rep(NA, p)
  #Obtain estimates for each of E(beta_i):
  for(i in 1:p) {
    #Calculate the sum of h(x)*weight
    sum_weight_beta <- t(beta_sampled[, i]) %*% weight
    
    #Estimate E(beta).
    E_beta_hat[i] <- sum_weight_beta / sum_weight
  }
  
  return(E_beta_hat)
}



###########################################
# Testing performance on repeated sampling
###########################################
estimates <- matrix(nrow = k, ncol = p)
for(i in 1:k) {
  YX <- generate_sample(beta)
  estimates[i, ] <- IS_estimates(YX[, 1], YX[, -1])
}
estimates

par(mfrow = c(1, p))
for(i in 1:p) {
  hist(estimates[, i], xlim = c(0.92*min(estimates[, i], beta[i]), 1.08*max(estimates[, i], beta[i])))
  abline(v = beta[i], col = "orange")
  abline(v = mean(estimates[, i]), col = "blue")
}



