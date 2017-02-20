
#Estimating P(X > a) where X ~ N(0, 1)
a <- 3
N <- 20000
m <- 300

estimate <- matrix(NA, nrow = m, ncol = 6)

#Standard MC approach:
x <- matrix(rnorm(N*m), nrow = m)
estimate[, 1] <- apply(x, 1, function(x) { sum(x > a) })/N

#IS with proposal distribution g(x) ~ t(df = 1)
x <- matrix(rt(N*m, 1), nrow = m)
estimate[, 2] <- apply(x, 1, function(x) { 
  sum(dnorm(x[x > a])/dt(x[x > a], 1)) })/N

#IS with proposal distribution g(x) ~ cauchy
x <- matrix(rcauchy(N*m), nrow = m)
estimate[, 3] <- apply(x, 1, function(x) { 
  sum(dnorm(x[x > a])/dcauchy(x[x > a])) })/N

#IS with proposal distribution g(x) ~ n(a, 1)
x <- matrix(rnorm(N*m, ifelse(a > 0, a, 0), 1), nrow = m)
estimate[, 4] <- apply(x, 1, function(x) { 
  sum(dnorm(x[x > a])/dnorm(x[x > a], ifelse(a > 0, a, 0), 1)) })/N

#IS with proposal dist g(x) ~ exp(1) + a (shifted exponential).
x <- matrix(rexp(N*m, 1) + a, nrow = m)
estimate[, 5] <- apply(x, 1, function(x) { sum(dnorm(x)/dexp(x - a)) })/N

#Change of variables: y = 1/x. Compute integral using unif(0, 1/a)
if(a == 0) { a <- 10^(-2)}
x <- matrix(runif(N*m, min(0, 1/a), max(0, 1/a)), nrow = m)
estimate[, 6] <- apply(x, 1, function(x) { abs(1/a)*sum(dnorm(1/x)/x^2) })/N
if(a < 0) { estimate[, 6] <- 1 - estimate[, 5] }


par(mfrow = c(2, 3))
xmin <- min(estimate)
xmax <- max(estimate)
true_val <- 1-pnorm(a)
apply(estimate, 2, function(x) { hist(x, xlim = c(xmin, xmax)); 
  abline(v = true_val, col = "orange") })


