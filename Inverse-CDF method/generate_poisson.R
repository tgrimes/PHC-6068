#Inverse-cdf for Poisson
rpoisson <- function(n = 1, lambda) {
  x <- rep(0, n)
  for(k in 1:n) {
    u <- runif(1)
    i <- floor(lambda + 0.33 - 0.02/lambda) #Approximate median.
    #Divide the search in half.
    if(u < ppois(i, lambda)) {
      #Search each intervals below the median.
      while(u < ppois(i - 1, lambda)) {
        i <- i - 1
      }
      x[k] <- i
    } else {
      #Search the intervals above the median.
      while(u > ppois(i, lambda)) {
        i <- i + 1
      }
      x[k] <- i
    }
  }
  return(x)
}

n <- 10000
lambda <- 3.7
y <- rpoisson(n, lambda)
hist(y, breaks = (min(y):(max(y) + 1) - 0.5),
     main = paste("Inverse-CDF method\nLambda =", lambda))
lines(min(y):max(y), n*dpois(min(y):max(y), lambda), 
      col = "red")
