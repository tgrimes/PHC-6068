library(rgl)
library(MASS)
library(scatterplot3d)
source("MCMC_Tools.R")

num_procedures <- 4
procedure_name <- matrix(0, nrow = num_procedures, ncol = 1)
cpu_times <- matrix(0, nrow = num_procedures, ncol = 3)
min_coordinates<- matrix(0, nrow = num_procedures, ncol = 2)
min_values <- matrix(0, nrow = num_procedures, ncol = 1)
X <- vector("list", num_procedures)

sim <- 30000

get.results <- function(digits = 4) {
  tab <- matrix(c(procedure_name, format(cpu_times, digits = 4), 
                  format(c(min_coordinates, min_values), digits = digits)), 
                nrow = num_procedures, byrow = F, 
                dimnames = list(NULL, c("Procedure", 
                                        "User", "System", "Total",
                                        "X", "Y", 
                                        "U(X, Y)")))
  return(as.data.frame(tab))
}
print.all <- function(digits = 4) {
  print(get.results(digits))
  print("Global minimum at (-+1.0445, -1.0084) with U = -8.12465")
}

#Function of interest.
# - point: point to evaluate. For a set of n points, an n by p matrix is expected.
# return: value of U at each point. If point is out of domain [-1.1, 1.1]^2, 
#         NA is returned.
U <- function(point) {
  n <- nrow(point) #Number of points.
  if(is.null(n)) {
    #If only 1 point, set up as a matrix.
    point <- matrix(point, nrow = 1, ncol = 2)
    n <- 1
  }
  x <- point[, 1]
  y <- point[, 2]
  sol <- -(x*sin(20*y) + y*sin(20*x))^2*cosh(sin(10*x)*x) -
            (x*cos(10*y) - y*sin(10*x))^2*cosh(cos(20*y)*y) 
  if(any(c(point < -1.1, point > 1.1))) {
    #warning("in U: point must be in [-1.1, 1.1]^2. Returning NA for these values.")
    if(n == 1) {
      index <- 1
    } else {
      index <- which(apply(point, 1, function(p) { any(c(p < -1.1, p > 1.1)) }))
    }
    for(i in index) {
      sol[i] <- NA
    }
  }
  return(sol)
}

#Graph of the function.
x <- seq(-1.1, 1.1, length.out = 300)
y <- seq(-1.1, 1.1, length.out = 300)
mesh <- expand.grid(x, y)
vals <- U(mesh)
plot3d(mesh[, 1], mesh[, 2], vals, size = 0.01, alpha = 0.15)

start_point <- c(0.5, 0.5)
##############################
# Attempt 1: MH algorithm
##############################

#Random walk proposal.
proposal <- proposal.random.walk(sigma = 0.001)

#Run random-walk MH and record CPU time.
procedure_name[1] <- "Random-walk MH"
cpu_times[1, ] <- as.numeric(system.time(
  X[[1]] <- MH_alg(U, proposal$q, proposal$rq, start_point, sim)
))[1:3]
vals <- U(X[[1]])
minimum <- which(vals == min(vals))[1]
min_coordinates[1, ] <- X[[1]][minimum, ]
min_values[1] <- vals[minimum]


#Independence proposal.
proposal <- proposal.independence.normal(mu = c(0, 0), sigma = 1)

#Run independence MH and record CPU time.
procedure_name[2] <- "Independence MH"
cpu_times[2, ] <- as.numeric(system.time(
  X[[2]] <- MH_alg(U, proposal$q, proposal$rq, start_point, sim)
))[1:3]
vals <- U(X[[2]])
minimum <- which(vals == min(vals))[1]
min_coordinates[2, ] <- X[[2]][minimum, ]
min_values[2] <- vals[minimum]




###########################################
# Attempt 2: Simulated annealing algorithm
###########################################
num_temp <- 30

#Random walk proposal.
proposal <- proposal.random.walk(sigma = 0.01)

#Run simmulated annealing with random-walk MH and log cooling.
t <- rep(10, num_temp)
for(i in 2:num_temp) { t[i:num_temp] <- log(t[i:num_temp] + 1)}
procedure_name[3] <- "annealing log"
cpu_times[3, ] <- as.numeric(system.time(
  X_temp <- annealing_alg(U, start_point, t, sigma = 1, n = sim)
))[1:3]
X[[3]] <- X_temp$X[[num_temp]]
vals <- U(X[[3]])
minimum <- which(vals == min(vals))[1]
min_coordinates[3, ] <- X[[3]][minimum, ]
min_values[3] <- vals[minimum]


#Run simmulated annealing with random-walk MH and square root cooling.
procedure_name[4] <- "annealing sqrt"
cpu_times[4, ] <- as.numeric(system.time(
  X_temp <- annealing_alg(U, start_point, 10^(1/2^(0:(num_temp - 1))) - 1,
                          sigma = 1, n = sim)
))[1:3]
X[[4]] <- X_temp$X[[num_temp]]
vals <- U(X[[4]])
minimum <- which(vals == min(vals))[1]
min_coordinates[4, ] <- X[[4]][minimum, ]
min_values[4] <- vals[minimum]




#Plotting the results:
plot_xy <- function(X, title) {
  plot(X[, 1], X[, 2], type = "l", lwd = 0.05, ylim = c(-1.2, 1.2), xlim = c(-1.2, 1.2),
       main = title, xlab = "x", ylab = "y")
  points(X[, 1], X[, 2], pch = 20, lwd = 0.1, cex = 0.1)
}
plot_x <- function(X) { 
  plot(1:length(X[, 1]), X[, 1], type = "l", lwd = 0.4, ylim = c(-1.2, 1.2),
       xlab = "iteration", ylab = "x") 
}


par(mfrow = c(1, 1))
plot_xy(X[[1]], procedure_name[[1]])

levels <- 10 #Number of contour lines.
par(mfrow = c(2, 3))
vals_U <- U(mesh)
min_U <- min(vals_U); max_U <- max(vals_U)
contour(x, y, matrix(vals_U, nrow = 300, ncol = 300),
        levels = min_U + (max_U - min_U)*seq(0, 1, 1/levels), 
        main = "Contours of U")
plot_xy(X[[1]], procedure_name[[1]])
plot_xy(X[[2]], procedure_name[[2]])
vals_psi <- exp(-U(mesh)/0.1)
min_psi <- min(vals_psi); max_psi <- max(vals_psi)
quantiles <- quantile(vals_psi, 1 - c(0.001, 0.99))
contour(x, y, matrix(vals_psi, nrow = 300, ncol = 300),
        levels = min_psi + (max_psi - min_psi)*seq(0, 1, 1/levels),
        main = "Contours of exp(-U(x)/t)")
plot_xy(X[[3]], procedure_name[[3]])
plot_xy(X[[4]], procedure_name[[4]])


print.all(5)

starting_point <- c(0.01, 0.01)
X_temp <- annealing_alg(U, starting_point, t, sigma = 1, n = 20000)
par(mfrow = c(3, 3))
index <- c(1:5, 7, 15, 20, 35)
for(i in 1:9) {
  j <- index[i]
  plot_xy(X_temp$X[[j]], paste("temp =", (X_temp$temps)[j]))
}
X[[3]] <- X_temp$X[[num_temp]]
vals <- U(X[[3]])
minimum <- which(vals == min(vals))[1]
min_coordinates[3, ] <- X[[3]][minimum, ]
min_values[3] <- vals[minimum]
print.all(5)

