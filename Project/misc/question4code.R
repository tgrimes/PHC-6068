set.seed(10)

# #####################################################
#
# Auxillary functions.
#
# #####################################################
#Generate n cities from a height by width square mile region.
generate_cities <- function(n = 100, height = 100, width = 100) {
  cities <- matrix(0, nrow = n, ncol = 2)
  cities[, 1] <- runif(n, 0, width)
  cities[, 2] <- runif(n, 0, height)
  
  return(cities)
}

#Swap two elements in a vector.
swap <- function(x, move) {
  move <- sort(move)
  x[seq(move[1], move[2])] <- x[seq(move[2], move[1])]
  x
}

#Distance measure between two cities.
distance <- function(city_a, city_b) {
  d <- city_b - city_a
  dist <- t(d) %*% d
  sqrt(dist)
}

#Create a scatterplot representing a tour.
tour_plot <- function(tour, cities, POINTS = TRUE, ...) {
  len <- round(tour_length(tour, cities), digits = 2)
  plot(c(cities[tour, 1], cities[tour[1], 1]), 
       c(cities[tour, 2], cities[tour[1], 2]),
       type = "l", main = paste("tour length =", len, "miles"),
       xlab = "", ylab = "", ...)
  if(POINTS) {
    color <- hsv(seq(0, 1, 1/n), 1, 0.95)
    points(cities[tour, 1], cities[tour, 2], pch = 16, col = color)
  }
}

#Compute the length of a tour.
tour_length <- function(tour, cities) {
  #Number of cities
  n <- length(tour)
  
  #Compute length of tour by summing the distances between each city.
  d <- 0
  for(i in 1:(n-1)) {
    d <- d + distance(cities[tour[i], ], cities[tour[i+1], ]) #Distance between cities i and i+1.
  }
  d <- d + distance(cities[tour[n], ], cities[tour[1], ]) #Loop end of tour back to beginning.
  
  d
}

#Compute the change in length of a tour after swapping cities i and j.
tour_change <- function(tour, cities, move) {
  #return(tour_length(swap(tour, move), cities) - tour_length(tour, cities))
  n <- nrow(cities)
  
  i <- move[1]
  j <- move[2]
  im1 <- ifelse(i - 1 < 1, n, i - 1)
  jp1 <- ifelse(j + 1 > n, 1, j + 1)
  
  new_distance <- distance(cities[tour[im1], ], cities[tour[j], ]) +
    distance(cities[tour[i], ], cities[tour[jp1], ])
  
  old_distance <- distance(cities[tour[im1], ], cities[tour[i], ]) +
    distance(cities[tour[j], ], cities[tour[jp1], ])

  new_distance - old_distance
}



# #####################################################
#
# Simulated annealing.
#
# #####################################################
#Simulated Annealing algorithm.
SA <- function(cities, tour = NULL, initial_t = 10, ITER = 10^5, HISTORY = TRUE) {
  n <- nrow(cities) #Number of cities
  
  #Initialize a tour.
  if(is.null(tour)) {
    tour <- sample(1:n, n)
  }
  
  #Initialize the temperature
  t <- initial_t
  
  #Keep track of all the tours (for plotting purposes later).
  if(HISTORY) {
    tour_history <- matrix(0, nrow = ITER/2, ncol = n)
    tour_history[1, ] <- tour
    changes <- 2
  }
  
  #Begin SA:
  for(iter in 1:ITER) {
    #Sample a new tour.
    move <- sort(sample(1:n, 2)) #Lin move. Choose two points to swap
    new_tour <- swap(tour, move)
    
    #Compute the change in distance of the tour length.
    delta_dist <- tour_change(tour, cities, move)
    
    #Determine if new tour should be accepted.
    if(delta_dist < 0) {
      tour <- new_tour
      if(HISTORY) {
        tour_history[changes, ] <- tour
        changes <- changes + 1
      }
    } else {
      u <- runif(1, 0, 1)
      if(exp(-delta_dist/t) > u) {
        tour <- new_tour
        if(HISTORY) {
          tour_history[changes, ] <- tour
          changes <- changes + 1
        }
      }
    }
    
    #Cool the temperature (square root cooling rate).
    t <- initial_t/sqrt(iter)
  }
  
  if(HISTORY) {
    tour_history = tour_history[1:(changes-1), ]
  } else {
    tour_history = NULL
  }
  
  return(list(shortest_tour = tour_length(tour, cities),
              tour_history = tour_history))
}


# #####################################################
#
# Stochastic approximation annealing.
#
# #####################################################
#Find which subset t belongs to in the sampling space partition.
get_E_index<-function(t, E){
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

#Simulated Annealing algorithm.
SAA <- function(cities, tour = NULL, initial_t = 10, ITER = 10^5, HISTORY = TRUE) {
  
  n <- nrow(cities) #Number of cities
  
  # Define initial values for SAMC
  E <- seq(750, 1000, length.out = 100)
  m <- length(E)
  theta <- rep(0, m + 1)
  PI <- 1/(m + 1) #(m + 1):1/(sum(1:(m + 1)))
  t0 <- 1000
  gamma <- t0/pmax(t0, 1:ITER)
  
  #Initialize a tour.
  if(is.null(tour)) {
    tour <- sample(1:n, n)
  }
  
  #Compute initial tour length.
  dist <- tour_length(tour, cities)
  
  #Find E index of current tour length.
  index <- get_E_index(dist, E)
  
  #Initialize the temperature
  t <- initial_t
  
  #Keep track of all the tours (for plotting purposes later).
  if(HISTORY) {
    tour_history <- matrix(0, nrow = ITER/2, ncol = n)
    tour_history[1, ] <- tour
    changes <- 2
  }
  
  #Begin SAA:
  for(iter in 1:ITER) {
    #Sample a new tour.
    move <- sort(sample(1:n, 2)) #Lin move. Choose two points to swap
    new_tour <- swap(tour, move)
    
    #Compute the change in distance of the tour length.
    delta_dist <- tour_change(tour, cities, move)

    #Determine if new tour should be accepted.
    if(delta_dist < 0) { 
      tour <- new_tour
      dist <- dist + delta_dist
      index <- get_E_index(dist, E)
      if(HISTORY) {
        tour_history[changes, ] <- tour
        changes <- changes + 1
      }
    } else {
      u <- runif(1, 0, 1) 
      new_index <- get_E_index(dist + delta_dist, E)
      if(exp(-delta_dist/t + theta[index] - theta[new_index]) > u){
        tour <- new_tour
        dist <- dist + delta_dist
        index <- new_index
        if(HISTORY) {
          tour_history[changes, ] <- tour
          changes <- changes + 1
        }
      } 
    }
    
    #Update theta.
    indicator <- rep(0, m + 1)
    indicator[index] <- 1
    theta <- theta + gamma[iter]*(indicator - PI)
    
    #Cool the temperature (square root cooling rate).
    t <- initial_t/sqrt(iter)
  }
  
  
  if(HISTORY) {
    tour_history = tour_history[1:(changes-1), ]
  } else {
    tour_history = NULL
  }
  
  return(list(shortest_tour = tour_length(tour, cities),
              tour_history = tour_history,
              theta = theta))
}



# #####################################################
#
# Simulations.
#
# #####################################################
#m simulations are performed; the shortest tour from SA and SAA will be saved.
m <- 100
shortest_tours <- matrix(0, nrow = m, ncol = 2)

##########################################
# Perform one simulation and store graphs.
##########################################
n <- 100 
cities <- generate_cities(n)

#Run SA and create graphs.
SA_results <- SA(cities, tour = 1:n)

png("images/SA_TSP_history.png", 1800, 1800, res = 240)
par(mfrow = c(3, 3), mar = rep(2, 4))
for(index in ceiling(seq(1, nrow(SA_results$tour_history), length.out = 9))) {
  tour_plot(SA_results$tour_history[index, ], cities, axes = FALSE)
  print(SA_results$shortest_tour)
}
dev.off()

#Run SAA and create graphs.
SAA_results <- SAA(cities, tour = 1:n)

png("images/SAA_TSP_history.png", 1800, 1800, res = 240)
par(mfrow = c(3, 3), mar = rep(2, 4))
for(index in ceiling(seq(1, nrow(SAA_results$tour_history), length.out = 9))) {
  tour_plot(SAA_results$tour_history[index, ], cities, axes = FALSE)
  print(tour_length(SAA_results$tour_history[index, ], cities))
}
dev.off()

shortest_tours[1, ] <- c(SA_results$shortest_tour, SAA_results$shortest_tour)

###################################
# Perform m - 1 more simulations.
###################################
for(i in 2:m) {
  shortest_tours[i, 1] <- (SA(generate_cities(n), tour = 1:n, 
                              HISTORY = FALSE))$shortest_tour
  shortest_tours[i, 2] <- (SAA(generate_cities(n), tour = 1:n, 
                               HISTORY = FALSE))$shortest_tour
}

write.csv(shortest_tours, "shortest_tours_10^5_iterations")
library("ggplot2")
colnames(shortest_tours) <- c("SA", "SAA")

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 1)
png("images/SA_vs_SAA.png", 400, 400, res = 100)
plot(shortest_tours[, 1], shortest_tours[, 2], pch = 16, 
     main = "Shortest paths from in 100 simulations",
     xlab = "SA", ylab = "SAA")
abline(v = mean(shortest_tours[, 1]), col = "red")
abline(h = mean(shortest_tours[, 2]), col = "red")
dev.off()
