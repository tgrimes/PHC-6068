#MCMC Tools

############################################
#
# Algorithms.
#
############################################
#Metropolis-Hastings algorithm
# - pi: function of interest.
# - q: density for proposal distribution
# - rq: function to generate random variables from proposal.
# - starting_value: initial point in sample space.
# - m: number of samples to generate.
# return: a sample of size n from pi.
MH_alg <- function(pi, q, rq, starting_value, m = 10^4) {
  #M-H ratio:
  ratio <- function(X, Y) {
    #Check that both numerator and denominator are numeric, 
    #and check for divide by zero.
    num <- pi(Y)*q(Y, X)
    den <- pi(X)*q(X, Y)
    if(any(c(is.na(den), is.na(num), den == 0))) {
      return(0)
    }
    num/den
  }
  
  X <- matrix(0, nrow = m, ncol = length(starting_value))
  X[1, ] <- starting_value
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

#Simulated annealing algorithm
# - H: function to minimize.
# - starting_value: initial point in sample space.
# - temperatures: vector of temperatures.
# - sigma: standard deviation for sampler, if none provided.
# - q: density for proposal distribution. DEFAULT is random-walk.
# - rq: function to generate random variables from proposal.
# - n: number of iterations at each temperature (can be a vector).
# return: the sampled points from final iteration.
annealing_alg <- function(H, starting_value, temperatures, sigma = NULL, 
                          q = NULL, rq = NULL, n = 10^4) {
  num_temps <- length(temperatures)
  dim <- length(starting_value)
  #Temperatures must be in decreasing order.
  temperatures <- temperatures[order(temperatures, decreasing = TRUE)]
  if(temperatures[num_temps] > 0.1) {
    warning("in annealing_alg: temperatures should approach 0.")
  }
  if(any(temperatures <= 0)) {
    warning("in annealing_alg: temperatures must be positive.")
  }
  if(is.null(q)) {
    if(is.null(sigma)) {
      stop("sigma required when q not provided.")
      # min_range <- 100 #Largest range to consider.
      # for(d in domain) {
      #   domain_range <- diff(d)
      #   min_range <- ifelse(domain_range < min_range, 
      #                       domain_range, min_range)
      # }
      # sigma <- 0.1*min_range*temperatures/max(temperatures)
    } else if(length(sigma) == 1) {
      sigma <- 0.1*sigma*temperatures/max(temperatures)
    }
    
    proposal <- proposal.random.walk(sigma)
    q <- vector("list", num_temps)
    rq <- vector("list", num_temps)
    for(i in 1:num_temps) {
      q[[i]] <- proposal[[i]]$q
      rq[[i]] <- proposal[[i]]$rq
    }
  } else if(length(q) == 1) {
    temp_q <- q
    temp_rq <- rq
    q <- vector("list", num_temps)
    rq <- vector("list", num_temps)
    for(i in 1:num_temps) {
      q[[i]] <- temp
      rq[[i]] <- temp
    }
  }
  
  if(length(n) == 1) {
    n <- rep(n, num_temps) #Ensure n is a vector.
  }
  #Boltzmann distribution.
  # point: the point to evaluate at (can be a vector).
  # U: the energy function to use.
  #return: the function exp(-H(x)/t)
  boltzmann <- function(t) {
    f <- function(x) { exp(-H(x)/t) }
    return(f)
  }
  
  X <- vector("list", num_temps)
  last_value <- starting_value
  for(i in 1:num_temps) {
    m <- n[i]
    X[[i]] <- matrix(0, nrow = m, ncol = dim)
    X[[i]] <- MH_alg(boltzmann(t[i]), q[[i]], rq[[i]], last_value, m)
    last_value <- X[[i]][m, ]
  }
  return(list(temps = temperatures, n = n, X = X))
}


#Simulated annealing algorithm
# - H: function to minimize.
# - q: density for proposal distribution
# - rq: function to generate random variables from proposal.
# - starting_value: initial point in sample space.
# - temperatures: vector of temperatures.
# - n: number of iterations at each temperature (can be a vector).
# return: the sampled points from final iteration.
# annealing_alg <- function(H, q, rq, starting_value, temperatures, n = 10^4) {
#   num_temps <- length(temperatures)


############################################
#
# Proposal Distributions.
#
############################################
#Random walk proposal:
#We will use an random-walk proposal T(X, Y) = q(Y) with Y|X ~ N(X, sigma)
proposal.random.walk <- function(sigma) {
  if(length(sigma) > 1) {
    k <- length(sigma)
    proposal <- vector("list", k)
    for(i in 1:k) {
      proposal[[i]] <- proposal.random.walk(sigma[i])
    }
    return(proposal)
  }
  
  sigma <- diag(2)*sigma
  sigma2_inv <- solve(sigma^2)
  #Density of the proposal (up to a proportion).
  q <- function(X, Y) {
    Y <- Y 
    mu <- X
    exp(-0.5*(t(Y - mu) %*% sigma2_inv %*% (Y - mu)))
  }
  #Generator for the proposal distribution.
  rq <- function(n = 1, X) {
    mu <- X
    mvrnorm(n, mu, sigma)
  }
  
  return(list(q = q, rq = rq, params = list(sigma = sigma), 
              name = "random-walk"))
}


#Independence Sampler
#Proposal distribution T(X -> Y):
#We will use an independence sampler T(X, Y) = q(Y), with Y|X ~ N(mean, sigma)
proposal.independence.normal <- function(mu, sigma) {
  if(!is.matrix(sigma)) {
    sigma <- 1.5*diag(2)*sigma #If covariance matrix is not given, create it.
  }
  sigma2_inv <- solve(sigma^2)
  #Density of the proposal (up to a proportion).
  q <- function(X, Y) {
    exp(-0.5*(t(Y - mu) %*% sigma2_inv %*% (Y - mu)))
  }
  #Generator for the proposal distribution.
  rq <- function(n = 1, X) {
    mvrnorm(n, mu, sigma)
  }
  
  return(list(q = q, rq = rq, params = list(mu = mu, sigma = sigma),
              name = "independence"))
}




