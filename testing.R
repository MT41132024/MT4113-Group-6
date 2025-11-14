### First Draft of functions
library(numDeriv)

# Univariate Newton Method
UVN <- function(start, f, data, epsilon){
  beta <- start
  delta <- 10^3
  while (abs(delta/beta) > epsilon){
    g <- grad(f, beta, data = data)
    h <- hessian(f, beta, data = data)[1,1]
    
    delta <- -g/h
    beta <- beta + delta
  }
  return(beta)
}

# Multivariate Newton Method
MVN <- function(start, f, data, epsilon){
  beta <- start
  delta <- 10^3
  # For each step, calculate gradient, hessian, solve linear system for delta
  while (min(abs(delta/beta)) > epsilon){ # Take minimum abs. change in each 
    g <- grad(f, beta, data = data)       # par, not total step size
    H <- hessian(f, beta, data = data)
    
    delta <- solve(H, -g)
    
    beta <- beta + delta
  }
  return(beta)
}

# Gauss-Newton Method (nls)
GN <- function(start, f, data, epsilon){
  theta <- start
  delta <- 10^3
  
  while (min(abs(delta/theta)) > epsilon) {
    j <- jacobian(function(x, data) data$y - f(data$x, x), theta, data = data)
    r <- data$y - f(data$x, theta)
    
    g <- t(j) %*% r
    H <- t(j) %*% j
    
    delta <- solve(H, -g)
    
    theta <- theta + delta
  }
  return(theta)
}
