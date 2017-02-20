set.seed(1000)
library(gridExtra)

# #####################################################
#
# Auxillary functions.
#
# #####################################################
#Generate data set for t-test.
generate_data <- function(n = 100, block = TRUE, weights = c(2, 2, 1)) {
  if(block) {
    n <- n/2
    W_top <- matrix(round(runif(n^2)*weights[1]), nrow = n, ncol = n)
    W_bottom <- matrix(round(runif(n^2)*weights[2]), nrow = n, ncol = n)
    epsilon <- matrix(runif(n^2)*weights[3], nrow = n, ncol = n)
    W <- rbind(cbind(W_top + t(W_top), epsilon),
               cbind(epsilon, W_bottom + t(W_bottom)))
    diag(W) <- 0
    return(W)
  } else {
    W <- matrix(0, nrow = n, ncol = n)
    for(i in 1:n) {
      #if(runif(1) > 0.1) {
        m <- rbinom(1, n, 0.05)
        latent_vars <- sample((1:n)[-i], m)
        vals <- runif(m)*weights[1]
        W[i, latent_vars] <- W[i, latent_vars] + vals
        W[latent_vars, i] <- W[i, latent_vars] + vals
      #}
    }
    return((W + t(W))/2)
  }
  # out_degree <- rnbinom(n, ceiling(0.1*n), mu = ceiling(0.1*n))
  # out_degree[which(out_degree < 1)] <- 1
  # out_degree[which(out_degree > (n - 1))] <- n - 1
  # graph <- get.adjacency(sample_degseq(out_degree, method = "simple.no.multiple"))
  # graph
  # return(graph)
  # plot(graph)
  
}

# W: 2n by 2n adjacency matrix
# A: vector of size n containing the indicies for one partition.
# a: include if cost of only one node is desired; assumed that a is in A.
external_cost <- function(W, A, a = NULL) {
  if(!is.null(a)) {
    return(sum(W[a, -A]))
  }
  c <- rep(0, nrow(W))
  c[A] <- apply(W[A, -A], 1, sum)
  c[-A] <- apply(W[-A, A], 1, sum)
  return(c)
}

# W: 2n by 2n adjacency matrix
# A: vector of size n containing the indicies for one partition.
# a: include if cost of only one node is desired; assumed that a is in A.
internal_cost <- function(W, A, a = NULL) {
  if(!is.null(a)) {
    return(sum(W[a, A]))
  }
  c <- rep(0, nrow(W))
  c[A] <- apply(W[A, A], 1, sum)
  c[-A] <- apply(W[-A, -A], 1, sum)
  return(c)
}

#Swap two nodes between partitions A and -A.
move <- function(W, A) {
  B <- (1:(2*length(A)))[-A]
  a <- sample(A, 1)
  b <- sample(B, 1)
  delta_t <- internal_cost(W, A, a) + internal_cost(W, B, b) -
    external_cost(W, A, a) - external_cost(W, B, b) + 2*W[a, b]
  A_new <- A
  A_new[which(A == a)] <- b
  return(list(A_new = A_new, delta_t = delta_t))
}

cut_size <- function(W, A) {
  if(length(unique(A)) != length(A)) {
    stop("A should contain unique elements.")
  }
  return(sum(external_cost(W, A)[A]))
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



## Approximate bisection
# returns a bisection of the graph that minimizes the cost using Kernighan/Lin Algorithm
# http://www.eecs.berkeley.edu/~demmel/cs267/lecture18/lecture18.html#link_4.2
# partition<-approximateBisection(W)
# W is symmetric matrix of size 2Nx2N made of non-negative values.
# partition is a list of two vectors of N indices.
# C.Ladroue
# #####################################################
#
# KL
#
# #####################################################
KL <- function(W, A = NULL){
  #minimumGain <- 1e-5 
  #minimum value for gain,
  # setting it to 0 might lead to infinite loop due to numerical inaccuracy
  
  N <- dim(W)[1] # number of elements
  m <- N/2
  
  # start off with a random partition
  if(is.null(A)) {
    A <- sample(1:N,N/2,replace=FALSE)
  }
  B <- (1:N)[-A]
  
  DA <- rowSums(W[A,B]) - rowSums(W[A,A]) + diag(W[A,A])
  DB <- rowSums(W[B,A]) - rowSums(W[B,B]) + diag(W[B,B])
  unmarkedA<-1:m
  unmarkedB<-1:m
  markedA<-rep(0,m)
  markedB<-rep(0,m)
  gains<-rep(0,m)
  for(k in 1:m){
    
    dimension<-m+1-k
    fasterGain<-matrix(DA[unmarkedA],nrow=dimension,ncol=dimension,byrow=FALSE)+
      matrix(DB[unmarkedB],nrow=dimension,ncol=dimension,byrow=TRUE)-
      2*W[A[unmarkedA],B[unmarkedB]]
    # mark the best pair
    best<-arrayInd(which.max(fasterGain),.dim=c(dimension,dimension))
    besti<-unmarkedA[best[1]]
    bestj<-unmarkedB[best[2]]
    bestGain<-fasterGain[best]
    markedA[k]<-unmarkedA[best[1]]
    markedB[k]<-unmarkedB[best[2]]
    unmarkedA<-unmarkedA[-best[1]]
    unmarkedB<-unmarkedB[-best[2]]
    # record gain
    gains[k]<-bestGain
    
    # update D for unmarked indices 
    DA[unmarkedA]<-DA[unmarkedA]+2*W[A[unmarkedA],A[besti]]-2*W[A[unmarkedA],B[bestj]]
    DB[unmarkedB]<-DB[unmarkedB]+2*W[B[unmarkedB],B[bestj]]-2*W[B[unmarkedB],A[besti]]
  }
  gains<-cumsum(gains)
  bestPartition<-which.max(gains)
  maxGain<-gains[bestPartition]
  if(maxGain>0){ 
    # swap best pairs
    A1<-c(A[-markedA[1:bestPartition]],B[markedB[1:bestPartition]])
    B1<-c(B[-markedB[1:bestPartition]],A[markedA[1:bestPartition]])
    A<-A1
    B<-B1
  }
  
  list(A = A, B = B, max_gain = maxGain)
}


# #####################################################
#
# Simulated annealing.
#
# #####################################################
#Simulated Annealing algorithm.
SA <- function(W, A = NULL, initial_t = 50, MAX_ITER = 10^5) {
  if(nrow(W) %% 2 != 0) {
    stop("W should have an even number of nodes.")
  }
  n <- nrow(W)/2 
  if(is.null(A)) {
    A <- sample(1:(2*n), n)
  }
  
  #Compute cut size of initial partition.
  t <- cut_size(W, A)
  temp <- initial_t
  
  #Maintain history of minimum cut sizes.
  t_history <- rep(0, MAX_ITER/10)
  min_t_history <- rep(0, MAX_ITER/10)
  t_history[1] <- t
  min_t_history[1] <- t
  min_t <- t
  min_A <- A
  
  #Begin SA:
  for(i in 1:MAX_ITER) {
    changes <- move(W, A)
    A_new <- changes$A_new
    delta_t <- changes$delta_t
    t_new <- t + delta_t
    
    #Determine if new tour should be accepted.
    if(delta_t < 0) {
      A <- A_new
      t <- t_new
    } else if(exp(-delta_t/temp) > runif(1, 0, 1)) {
      A <- A_new
      t <- t_new
    }
    
    #Cool the temperature (square root cooling rate).
    temp <- initial_t/sqrt(i)
    
    #Update history of cut sizes.
    if(t_new < min_t) {
      min_A <- A      
      min_t <- t_new
    }
    if(i %% 10 == 0) {
      t_history[i/10] <- t
      min_t_history[i/10] <- min_t
    }
  }
  
  return(list(A = min_A, t_history = t_history, 
              min_t_history = min_t_history, 
              iterations = MAX_ITER))
}


# #####################################################
#
# Stochastic approximation annealing.
#
# #####################################################
SAA <- function(W, A = NULL, initial_t = 60, MAX_ITER = 10^4,
                E = seq(10000, 12000, length.out = 100),
                PI = NULL, t0 = 1000) {
  
  if(nrow(W) %% 2 != 0) {
    stop("W should have an even number of nodes.")
  }
  n <- nrow(W)/2 
  if(is.null(A)) {
    A <- sample(1:(2*n), n)
  }
  m <- length(E)
  
  # Define initial values for SAA
  #PI <- (m + 1):1/((m + 2)*(m + 1)/2)
  PI <- exp(-0.05*(1:(m + 1) - 1)); PI <- PI / sum(PI)
  #PI <- rep(1/(m + 1), m + 1)
  gamma <- t0/(pmax(t0, 1:MAX_ITER))
  theta <- rep(0, m + 1)
  
  #Compute cut size of initial partition.
  t <- cut_size(W, A)
  
  #Find E index of current cut size.
  index <- get_partition_index(t, E)
  
  #Initialize the temperature
  temp <- initial_t
  
  #Maintain history of minimum cut sizes.
  t_history <- rep(0, MAX_ITER/10)
  min_t_history <- rep(0, MAX_ITER/10)
  t_history[1] <- t
  min_t_history[1] <- t
  min_t <- t
  min_A <- A
  
  #Begin SAA:
  for(i in 1:MAX_ITER) {
    #Sample a new tour.
    changes <- move(W, A)
    A_new <-changes$A_new
    delta_t <- changes$delta_t
    t_new <- t + delta_t
    
    #Determine if new tour should be accepted.
    if(delta_t < 0) { 
      A <- A_new
      t <- t_new
      index <- get_partition_index(t, E)
    } else {
      u <- runif(1, 0, 1) 
      new_index <- get_partition_index(t + delta_t, E)
      if(exp(-delta_t/temp + theta[index] - theta[new_index]) > u){
        A <- A_new
        t <- t_new
        index <- new_index
      } 
    }
    
    #Update theta.
    indicator <- rep(0, m + 1)
    indicator[index] <- 1
    theta <- theta + gamma[i]*(indicator - PI)
    
    #Cool the temperature (square root cooling rate).
    temp <- initial_t/sqrt(i)
    
    #Update history of cut sizes.
    if(t_new < min_t) {
      min_A <- A      
      min_t <- t_new
    }
    if(i %% 10 == 0) {
      t_history[i/10] <- t
      min_t_history[i/10] <- min_t
    }
  }
  
  return(list(A = min_A, t_history = t_history, 
              min_t_history = min_t_history, 
              iterations = MAX_ITER))
}

# #####################################################
#
# KL algorithm
#
# #####################################################
#Returns the right tail probability at t.
# KL <- function(W, A = NULL) {
#   n <- nrow(W)
#   if(n %% 2 != 0) {
#     stop("W should have an even number of verticies.")
#   }
#   if(is.null(A)) {
#     A <- 1:(n/2)
#   }
#   
#   unlocked <- 1:n
#   
#   E <- external_cost(W, A)
#   I <- internal_cost(W, A)
#   D <- E - I
# }



# #####################################################
#
# SAMC
#
# #####################################################
SAMC <- function(W, A = NULL, MAX_ITER = 10^5, E = seq(1100, 1400, length.out = 100),
                 PI = NULL, t0 = 5000) {
  if(nrow(W) %% 2 != 0) {
    stop("W should have an even number of nodes.")
  }
  n <- nrow(W)/2 
  if(is.null(A)) {
    A <- sample(1:(2*n), n)
  }
  m <- length(E)
  
  #Initialize variables for SAMC.
  #PI <- (m + 1):1/((m + 2)*(m + 1)/2)
  PI <- exp(-0.05*(1:(m + 1) - 1)); PI <- PI / sum(PI)
  #PI <- rep(1/(m + 1), m + 1)
  gamma <- t0/(pmax(t0, 1:MAX_ITER))
  theta <- rep(0, m + 1)
  
  #Find the sampling partition index of the current cut size.
  t <- cut_size(W, A)
  partition_index <- get_partition_index(t, E)
  
  #Maintain history of minimum cut sizes.
  t_history <- rep(0, MAX_ITER/10)
  min_t_history <- rep(0, MAX_ITER/10)
  t_history[1] <- t
  min_t_history[1] <- t
  min_t <- t
  min_A <- A
  
  for(i in 2:MAX_ITER) {
    #Generate a permutation from the given sample.
    changes <- move(W, A)
    A_new <-changes$A_new
    delta_t <- changes$delta_t
    
    #Find the partition index of this t-statistic.
    t_new <- t + delta_t
    partition_index_new <- get_partition_index(t_new, E)
    
    #Compute the MH ratio.
    r <- exp(-delta_t + theta[partition_index] - theta[partition_index_new])
    
    if(r > runif(1)) {
      A <- A_new
      t <- t_new
      partition_index <- partition_index_new
    } else {
      #Do nothing.
    }
    
    #Update theta (based on whether this cut size was accepted).
    indicator <- rep(0, m + 1)
    indicator[partition_index] <- 1
    theta <- theta + gamma[i]*(indicator - PI)
    
    #Update history of cut sizes.
    if(t_new < min_t) {
      min_A <- A      
      min_t <- t_new
    }
    if(i %% 10 == 0) {
      t_history[i/10] <- t
      min_t_history[i/10] <- min_t
    }
    
  }
  
  return(list(A = min_A, t_history = t_history, 
              min_t_history = min_t_history, 
              iterations = MAX_ITER, theta = theta, PI = PI, E = E))
}


# #####################################################
#
# Simulations
#
# #####################################################
# library(igraph)
# Y <- W
# Y[A, -A] <- 0
# Y[-A, A] <- 0
# Y[A, A] <- 1
# Y[-A, -A] <- 1
# membership <- rep(0, n)
# membership[A] <- 1
# #g1 <- graph_from_adjacency_matrix(W, weighted = TRUE, mode = "undirected")
# #E(g1)$color
# #group <- groups(list(membership = membership))
# #tkplot(g1)
# 
# 
# 
n <- c(100, 250, 500, 750, 1000, 2000, 5000, 10000)
cutsize <- matrix(0, nrow = 4, ncol = length(n))
cuttime <- matrix(0, nrow = 4, ncol = length(n))
for(i in 1:length(n)) {
  iter <- 10^5
  MAX_TIME <- 100
  W <- generate_data(n = n[i], block = FALSE)
  # ##############
  # KL
  # ##############
  if(n[i] > 2000) {
    cutsize[1, i] <- NA
  } else {
    start_time <- proc.time()
    result1 <- KL(W)
    timeA1 <- proc.time() - start_time
    while(result1$max_gain > 1e-4 & timeA1[[3]] < MAX_TIME) {
      result1 <- KL(W, result1$A)
      timeA1 <- proc.time() - start_time
    }
    A1 <- result1$A
    cutsize[1, i] <- cut_size(W, A1)
    cuttime[1, i] <- timeA1[[3]]
  }
  
  # ##############
  # SA
  # ##############
  #Initial temperature is cutA1*0.1/log(2); this enables a 50% chance
  # of accepting a 10% upward move, initially.
  start_time <- proc.time()
  result2 <- SA(W, initial_t = 100, MAX_ITER = iter)
  timeA2 <- proc.time() - start_time
  A2 <- result2$A
  cutsize[2, i] <- cut_size(W, A2)
  cuttime[2, i] <- timeA2[[3]]
  
  # ##############
  # SAA
  # ##############
  start_time <- proc.time()
  result3 <- SAA(W, E = seq(result2$t_history[[1000]]*0.8, 
                            result2$t_history[[1000]]*1.5, length.out = 100), 
                 t0 = 5000, initial_t = 100, MAX_ITER = iter)
  timeA3 <- proc.time() - start_time
  A3 <- result3$A
  cutsize[3, i] <- cut_size(W, A3)
  cuttime[3, i] <- timeA3[[3]]
  # ##############
  # SAMC
  # ##############
  start_time <- proc.time()
  result4 <- SAMC(W, E = seq(result2$t_history[[1000]]*0.8, 
                             result2$t_history[[1000]]*1.5, length.out = 100),
                  t0 = 5000, MAX_ITER = iter)
  timeA4 <- proc.time() - start_time
  A4 <- result4$A
  cutsize[4, i] <- cut_size(W, A4)
  cuttime[4, i] <- timeA4[[3]]
  
  
  # #####################################################
  #
  # Results
  #
  # #####################################################
  png(paste("images/graph_min_cut_n", n[i], "_iter", iter, ".png", sep = ""), 
      2000, 2000, res = 400)
  colors <- c("black", "orange", "blue", "gray")
  plot((1:(result2$iterations/10))*10, result2$min_t_history, 
       type = "l",
       xlim = c(0, iter),
       ylim = c(min(cutsize[2:4, i]), max(result2$min_t_history)),
       ylab = "Minimum cut size", xlab = "Iteration number",
       col = colors[1], lty = 1, xaxt = "n")
  axis(1, at = 0:2*iter/2)
  lines((1:(result3$iterations/10))*10, result3$min_t_history, 
        col = colors[2], lty = 2)
  lines((1:(result4$iterations/10))*10, result4$min_t_history, 
        col = colors[3], lty = 3)
  legend("topright", c("SA", "SAA", "SAMC"), cex = 0.8, col = colors, 
         lty = 1:4, lwd = 2, bty = "n")
  dev.off()
  
  png(paste("images/graph_cut_n", n[i], "_iter", iter, ".png", sep = ""), 
      2000, 2000, res = 400)
  plot((1:(result2$iterations/10))*10, result2$t_history, 
       type = "l",
       xlim = c(0, iter),
       ylim = c(min(cutsize[2:4, i]), max(result2$t_history)),
       ylab = "Cut size", xlab = "Iteration number",
       col = colors[1], lty = 1, xaxt = "n")
  axis(1, at = 0:2*iter/2)
  lines((1:(result3$iterations/10))*10, result3$t_history, 
        col = colors[2], lty = 2)
  lines((1:(result4$iterations/10))*10, result4$t_history, 
        col = colors[3], lty = 3)
  legend("topright", c("SA", "SAA", "SAMC"), cex = 0.8, col = colors, 
         lty = 1:4, lwd = 2, bty = "n")
  dev.off()
}

png(paste("images/graph_all_vals_iter", iter, ".png", sep = ""), 
    600, 140, res = 100)
vals <- round(cutsize)
colnames(vals) <- n
rownames(vals) <- c( "KL", "SA", "SAA", "SAMC")
grid.table(vals)
dev.off()

png(paste("images/graph_all_times_iter", iter, ".png", sep = ""), 
    2000, 2000, res = 400)
plot(n[n <= 2000], cuttime[1, n <= 2000], type = "l", col = colors[4], lty = 4,
     ylim = c(0, max(cuttime)), ylab = "Total CPU time (s)", xlab = "n",
     xlim = c(n[1], n[length(n)]))
lines(n, cuttime[2, ], col = colors[1], lty = 1)
lines(n, cuttime[3, ], col = colors[2], lty = 2)
lines(n, cuttime[4, ], col = colors[3], lty = 3)
legend("topright", c("KL", "SA", "SAA", "SAMC"), cex = 0.8, 
       col = colors[c(4, 1, 2, 3)], lty = c(4, 1, 2, 3), lwd = 2, bty = "n")
dev.off()


cutsize
cuttime
