#Gelfand and Smith 1990, JASA
set.seed(10000)

I <- 30 #i = 1, ..., I   Number of rats.
J <- 5  #j = 1, ..., J   Number of measurements.

X <- c(8, 15, 22, 29, 36)
xbar <- mean(X)

#Hyperparameters.
a <- 0      #alpha_c ~ N(a, b)
b <- sqrt(1000)    #alpha_c ~ N(a, b)
c <- 0.01   #sigma_a,b,c^(-2) ~ Gamma(c, d)
d <- 0.01

#Conditional distributions. For each, y is expected to be the entire data set,
# y : I by J matrix containing the data set.
# x : J by 1 vector of x values.
# alpha : I by 1 vector of alpha values.
# beta : I by 1 vector of beta values.
# alpha_c, beta_c, sigma_a, sigma_b, sigma_c : all scalars.
cond_alpha <- function(Y, alpha_c, sigma_a_n2, sigma_c_n2) {
  A <- J * sigma_c_n2 + sigma_a_n2
  B <- (Y %*% rep(1, J)) * sigma_c_n2 + alpha_c * sigma_a_n2
  rnorm(I, B/A, A^(-0.5))
}
cond_beta<- function(Y, alpha, beta_c, sigma_b_n2, sigma_c_n2) {
  A <- sum((X - xbar)^2) * sigma_c_n2 + sigma_b_n2
  B <- (Y %*% (X - xbar)) * sigma_c_n2 + beta_c * sigma_b_n2
  rnorm(I, B/A, A^(-0.5))
}
cond_alpha_c <- function(alpha, sigma_a_n2) {
  A <- I * sigma_a_n2 + 1/b^2
  B <- sum(alpha) * sigma_a_n2 + a / b^2
  rnorm(1, B/A, A^(-0.5))
}
cond_beta_c <- function(beta, sigma_b_n2) {
  A <- I * sigma_b_n2 + 1/b^2
  B <- sum(beta) * sigma_b_n2 + a / b^2
  rnorm(1, B/A, A^(-0.5))
}
cond_sigma_a_n2 <- function(alpha, alpha_c) {
  A <- 0.5 * sum((alpha - alpha_c)^2) + d
  rgamma(1, I/2 + c, rate = A)
}
cond_sigma_b_n2 <- function(beta, beta_c) {
  A <- 0.5 * sum((beta - beta_c)^2) + d
  rgamma(1, I/2 + c, rate = A)
}
cond_sigma_c_n2 <- function(Y, alpha, beta) {
  Z <- Y - alpha - (beta %*% t(X - xbar))
  A <- 0.5 * sum(Z^2) + d
  rgamma(1, I*J/2 + c, rate = A)
}


#Gibbs sampler.
Gibbs <- function(Y, BURN = 100, ITER = 1000, init_alpha = rep(0, I), 
                  init_beta = rep(1, I), init_alpha_c = 1, init_beta_c = 1, 
                  init_sigma_a = 1, init_sigma_b = 1, init_sigma_c = 1,
                  GRAPH_DIAGNOSTIC = TRUE) {

  #Store estimates of alpha, beta, and csigma.
  alpha <- matrix(0, nrow = BURN + ITER, ncol = I)
  beta <- matrix(0, nrow = BURN + ITER, ncol = I)
  sigma_c_n2 <- rep(0, I)
  
  #Initialize estimates.
  alpha[1, ] <- init_alpha
  beta[1, ] <- init_beta
  sigma_c_n2[1] <- init_sigma_c^(-2)
  alpha_c <- init_alpha_c
  beta_c <- init_beta_c
  sigma_a_n2 <- init_sigma_a^(-2)
  sigma_b_n2 <- init_sigma_b^(-2)

  for(i in 2:(BURN + ITER)) {
    alpha[i, ] <- cond_alpha(Y, alpha_c, sigma_a_n2, sigma_c_n2[i - 1])
    beta[i, ] <- cond_beta(Y, alpha[i, ], beta_c, sigma_b_n2, sigma_c_n2[i - 1])
    sigma_c_n2[i] <- cond_sigma_c_n2(Y, alpha[i, ], beta[i, ])
    sigma_a_n2 <- cond_sigma_a_n2(alpha[i, ], alpha_c)
    sigma_b_n2 <- cond_sigma_b_n2(beta[i, ], beta_c)
    alpha_c <- cond_alpha_c(alpha[i, ], sigma_a_n2)
    beta_c <- cond_beta_c(beta[i, ], sigma_b_n2)
  }
  
  if(GRAPH_DIAGNOSTIC) {
    iterations <- seq(1, BURN, ceiling(0.005*BURN))
    png("images/alpha_diagnostic.png", 400, 400, res = 100)
    plot(x = iterations, y = alpha[iterations, 1], type = "l",
         main = "Gibbs sampling of alpha",
         xlab = "iteration number", ylab = "alpha")
    for(i in 1:I) {
      lines(x = iterations, y = alpha[iterations, i])
    }
    dev.off()
    
    png("images/beta_diagnostic.png", 400, 400, res = 100)
    plot(x = iterations, y = beta[iterations, 1], type = "l",
         main = "Gibbs sampling of beta",
         xlab = "iteration number", ylab = "beta")
    for(i in 1:I) {
      lines(x = iterations, y = beta[iterations, i])
    }
    dev.off()
    
    png("images/sigma_diagnostic.png", 400, 400, res = 100)
    plot(x = iterations, y = sigma_c_n2[iterations]^(-1/2), type = "l",
         main = "Gibbs sampling of sigma",
         xlab = "iteration number", ylab = "sigma")
    dev.off()
  }
  alpha_hat <- apply(alpha[(BURN + 1):(BURN + ITER), ], 2, mean)
  beta_hat <- apply(beta[(BURN + 1):(BURN + ITER), ], 2, mean)
  sigma_c_hat <- mean(sigma_c_n2[(BURN + 1):(BURN + ITER)]^(-1/2))
  
  return(list(alpha_hat = alpha_hat,
              beta_hat = beta_hat,
              sigma_c_hat = sigma_c_hat))
}

# #####################################################
#
# Part A.
#
# #####################################################
library(DPpackage)
data(rats)
rats$rat <- as.factor(rats$rat)
Y <- matrix(rats$weight, nrow = 30, ncol = 5, byrow = TRUE)

#Run Gibbs.
results <- Gibbs(Y, BURN = 10000)
results

#Graph the resulting regression lines onto the rats data.
Gibbsframe <- data.frame(rat = as.factor(1:30),
                         alpha_hat = results$alpha_hat,
                         beta_hat = results$beta_hat)

png("images/gibbs_regression.png", 680, 420, res = 100)
ggplot() +
  geom_point(data = rats, aes(day, weight, color = rat), alpha = 0.6) +
  geom_abline(data = Gibbsframe, alpha = 0.6,
              aes(intercept = alpha_hat - beta_hat*xbar, slope = beta_hat, color = rat))
dev.off()
# #####################################################
#
# Part B.
#
# #####################################################
#Generate a data set from specified parameter values.
#Note, sigma_a, sigma_b, and sigma_c are requested, not sigma^(-2).
generate_data <- function(alpha_c = 0, beta_c = 0, sigma_a = 1,
                          sigma_b = 1, sigma_c = 1) {
  #Generate alphas and betas.
  alpha <- rnorm(I, alpha_c, sigma_a)
  beta <- rnorm(I, beta_c, sigma_b)
  sigma_c <- abs(rnorm(1, sigma_c, 0.2*sigma_c))
                 
  #Generate y's.
  Y <- matrix(rnorm(I*J, alpha + beta %*% t(X - xbar), 
                    sigma_c), nrow = I, ncol = J)
  
  return(list(Y = Y, 
              params = list(alpha = alpha, beta = beta, 
                            sigma_c= sigma_c)))
}

#Compare the Gibbs estimates to the true parameters and summarize the results.
summarize <- function(results, params, GRAPH = TRUE) {
  #Analyze the results.
  true <- c(params$alpha, params$beta, params$sigma_c)
  pred <- c(results$alpha_hat, results$beta_hat, results$sigma_c_hat)
  label <- c(rep("alpha", I), rep("beta", I), "sigma")
  
  #Make scatterplot of relative errors.
  g <- NULL
  if(GRAPH) {
    g <- ggplot(data.frame(residuals = (pred - true)/true, rat = c(rep(1:I, 2), I/2), label = label),
                aes(rat, residuals), environment = environment()) +
      geom_point(aes(color = label), alpha = 0.7) +
      geom_abline(slope = 0, intercept = 0) +
      ylab("relative error")
    
    g
  }
  
  residuals <- (pred - true)/true
  alpha_residuals <- residuals[1:I]
  beta_residuals <- residuals[(I + 1):(2*I)]
  sigma_residual <- residuals[2*I + 1]
  
  return(list(alpha_residuals = alpha_residuals, 
              beta_residuals = beta_residuals,
              sigma_residual = sigma_residual,
              g = g))
}

#Compare the Gibbs estimates to the true parameters and summarize the results.
sim <- 1000
alpha_residuals <- rep(0, I*sim)
beta_residuals <- rep(0, I*sim)
sigma_residuals <- rep(0, sim)
for(i in 1:sim) {
  data <- generate_data(240, 6, 14, 0.5, 6)
  results <- Gibbs(data$Y, ITER = 10000, GRAPH_DIAGNOSTIC = FALSE)
  s <- summarize(results, data$params, GRAPH = FALSE)
  alpha_residuals[(1 + I*(i - 1)):(I*i)] <- s$alpha_residuals
  beta_residuals[(1 + I*(i - 1)):(I*i)] <- s$beta_residuals
  sigma_residuals[i] <- s$sigma_residual
  #ggsave(paste("images/results_rats_sim_", i, ".png", sep = ""))
}

png("images/alpha_rel_errors.png", 400, 400, res = 100)
hist(alpha_residuals, xlim = c(-0.2, 0.2), xlab = c("relative error"),
     main = paste("Relative error of alpha estimates\n",
                  "from", sim, "simulated data sets"))
dev.off()

png("images/beta_rel_errors.png", 400, 400, res = 100)
hist(beta_residuals, xlim = c(-0.2, 0.2), xlab = c("relative error"),
     main = paste("Relative error of beta estimates\n",
                  "from", sim, "simulated data sets"))
dev.off()

png("images/sigma_rel_errors.png", 400, 400, res = 100)
hist(sigma_residuals, xlim = c(-0.2, 0.2), xlab = c("relative error"),
     main = paste("Relative error of sigma estimates\n",
                  "from", sim, "simulated data sets"))
dev.off()


