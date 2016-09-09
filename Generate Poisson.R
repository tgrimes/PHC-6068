#Inverse-cdf for Poisson
rpoisson <- function(n = 1, lambda) {
  x <- rep(0, n)
  for(k in 1:n) {
    u <- runif(1)
    i <- floor(lambda)
    if(u < ppois(i, lambda)) {
      while(u < ppois(i - 1, lambda)) {
        i <- i - 1
      }
      x[k] <- i
    } else {
      while(u > ppois(i + 1, lambda)) {
        i <- i + 1
      }
      x[k] <- i
    }
  }
  return(x)
}


hist(rpoisson(100, 50))
