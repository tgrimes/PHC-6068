# Yu, Liang, Chatterjee, and Ciampa, 2011, Biostatistics
set.seed(1000)
library(gridExtra)

#Generate data set for t-test.
generate_data <- function(n = 1000, mu = c(0, 0), sigma = c(1, 1)) {
  Y <- cbind(rnorm(n, mu[1], sigma[1]),
             rnorm(n, mu[2], sigma[2]))
  colnames(Y) <- c("X", "Y")
  
  return(Y)
}

#Calculate t-test statistic for a given data set Y.
t_val <- function(Y) {
  if(any(is.na(Y))) {
    error("Y contains NA values.")
  }
  t <- as.numeric(t.test(Y[, 2], Y[, 1], var.equal = TRUE)$statistic)
}

# #####################################################
#
# Randomization test
#
# #####################################################
#Returns the right tail probability at t.
randomization_p_value<-function(t, t_sample){
  pval <- sum(t_sample[t_sample >= t])/length(t_sample)
  
  return(pval)
}

randomization <- function(Y, MAX_ITER = 10^5, 
                        PERMUTATION_PERCENT = 0.05) {
  n <- nrow(Y)
  #Number of observations to swap:
  swap_size <- ceiling(n*PERMUTATION_PERCENT) 
  t_sample <- rep(0, MAX_ITER)
  t_sample[1] <- t_val(Y)
  
  for(i in 1:MAX_ITER) {
    swap_indicies <- sample(1:n, swap_size)
    Y_new <- Y
    Y_new[swap_indicies, ] <- Y_new[swap_indicies, c(2, 1)]
    t_sample[i] <- t_val(Y_new)
    Y <- Y_new
  }
  
  return(t_sample)
}

# #####################################################
#
# SAMC
#
# #####################################################
#Returns the right tail probability P(T >= t)
SAMC_p_value<-function(t, E, theta, PI){
  m <- length(theta)
  psi <- exp(theta)*PI #Assuming no empty subregions.
  normalize <- sum(psi) #Normalizing factor.
  p <- psi/normalize
  pval <- 0
  
  index <- get_partition_index(t, E)
  if(index == 1) {
    pval <- 1
  } else {
    pval <- (t - E[index - 1])/
      (E[index] - E[index - 1])*sum(p[1:index])
    pval <- pval + (E[index] - t)/
      (E[index] - E[index - 1])*sum(p[1:(index - 1)])
    pval <- 1 - pval 
  }
  
  return(pval)
}

#Find which subset t belongs to in the sampling space partition.
get_partition_index<-function(t, E){
  m <- length(E)
  #First check edge cases.
  if(t < E[1]) {
    return(1)
  } else if (t >= E[m]) {
    return(m + 1)
  } 
  
  #Binary search to find which subset contains t.
  binary_search <- function(t, indicies) {
    m <- length(indicies)
    if(m == 1) {
      if(t < E[indicies]) {
        return(indicies)
      } else {
        return(indicies + 1)
      }
    }
    i <- ceiling(m/2) + 1 #Start in the middle of the set of indicies.
    if(t >= E[indicies[i - 1]] && t < E[indicies[i]]) {
      return(indicies[i])
    }
    if(t < E[indicies[i - 1]]) {
      return(binary_search(t, indicies[1:(i - 1)]))
    } else {
      return(binary_search(t, indicies[i:m]))
    }
  }
  
  return(binary_search(t, 2:m))
}

SAMC <- function(Y, MAX_ITER = 10^5, E = seq(0, 7, length.out = 100),
                 PI = NULL, t0 = 5000, PERMUTATION_PERCENT = 0.05) {
  n <- nrow(Y)   
  m <- length(E)
  #Number of observations to swap:
  swap_size <- ceiling(n*PERMUTATION_PERCENT) 
  
  #Initialize variables for SAMC.
  PI <- rep(1/(m + 1), m + 1)
  #PI <- 1:(m + 1)/((m + 2)*(m + 1)/2)
  #PI <- exp(-0.05*((m + 1):1 - 1)); PI <- PI / sum(PI)
  gamma <- t0/(pmax(t0, 1:MAX_ITER))
  theta <- rep(0, m + 1)
  
  #Find the sampling partition index of the current t-statistic.
  t <- t_val(Y)
  partition_index <- get_partition_index(t, E)
  
  for(i in 2:MAX_ITER) {
    #Generate a permutation from the given sample.
    swap_indicies <- sample(1:n, swap_size)
    Y_new <- Y
    Y_new[swap_indicies, ] <- Y_new[swap_indicies, c(2, 1)]
    
    #Find the partition index of this t-statistic.
    t_new <- t_val(Y_new)
    partition_index_new <- get_partition_index(t, E)
    
    #Compute the MH ratio.
    r <- exp(theta[partition_index] - theta[partition_index_new])
    
    if(r > runif(1)) {
      t <- t_new
      Y <- Y_new
      partition_index <- partition_index_new
    } else {
      #Do nothing.
    }
    
    #Update theta (based on whether this t-statistic was accepted).
    indicator <- rep(0, m + 1)
    indicator[partition_index] <- 1
    theta <- theta + gamma[i]*(indicator - PI)
  }
  
  return(list(theta = theta, PI = PI, E = E))
}

# #####################################################
#
# Results
#
# #####################################################
p <- 1*10^(-(2:10)) #pvalues to consider
randomization_iter <- 10^6
SAMC_iter <- 10^6
estimates <- matrix(0, nrow = length(p), ncol = 5)
colnames(estimates) <- c("True p-value", "Randomization", "SAMC", 
                         "ARE% (Rand)", "ARE% (SAMC)")
#Repeat simulation for 20 different data sets generated under H0.
Y <- generate_data()
n <- nrow(Y)

#Randomization estimate from sample.
# t_sample <- randomization(Y, randomization_iter)

#SAMC estimate from sample.
results <- SAMC(Y, SAMC_iter)

#Compute p-values for different t-statistic values
for(j in 1:length(p)) {
  #Obtain t statistic for desired p-value.
  t <- qt(1 - p[j], 1998)
  
  #Randomization estimate.
  randomization_val <- randomization_p_value(t, t_sample)
  randomization_ARE <- abs((randomization_val - p[j])/p[j])*100
  
  #SAMC estimates:
  SAMC_val <- SAMC_p_value(t, results$E, results$theta, results$PI)
  SAMC_ARE <- abs((SAMC_val - p[j])/p[j])*100
  
  estimates[j, ] <- 
    c(p[j], randomization_val, SAMC_val, randomization_ARE, SAMC_ARE)
}

# png("images/p_val_rand_dist.png", 400, 400, res = 100) {
#   hist(t_sample, main = "Distribution of t-statistics\n 
#        by randomization procedure",
#        xlab = "t-statistic", freq = FALSE)
#   curve(pt(x, 1998), min(t_sample), max(t_sample))
# }

tab1 <- estimates[, c(1, 3, 5)]
tab1[, 3] <- round(tab1[, 3], 1)
tab1[, 1] <- signif(tab1[, 1], digits = 1)
tab1[, 2] <- signif(tab1[, 2], digits = 3)
colnames(tab1) <- c("True P-value", "SAMC", "ARE%")
png("images/p_val_ARE.png", 400, 400, res = 100)
grid.table(tab1, rep("", nrow(tab1)))
dev.off()