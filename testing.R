### First Draft of functions
library(numDeriv)

# Univariate Newton Method
UVN <- function(f, inits, data, minimum=TRUE, tol=1e-5){
  theta <- inits
  delta <- 10^3
  
  if (minimum == FALSE){
    f <- -f
  }
  
  while (abs(delta/theta) > tol){
    g <- grad(f, theta, data = data)
    h <- hessian(f, theta, data = data)[1,1]
    
    delta <- -g/h
    theta <- theta + delta
  }
  return(theta)
}

# Multivariate Newton Method
MVN <- function(f, inits, data, minimum=TRUE, tol=1e-5){
  theta <- inits
  delta <- 10^3
  
  if (minimum == FALSE){
    f <- -f
  }
  
  # For each step, calculate gradient, hessian, solve linear system for delta
  while (min(abs(delta/theta)) > tol){ # Take minimum abs. change in each 
    g <- grad(f, theta, data = data)       # par, not total step size
    H <- hessian(f, theta, data = data)
    
    delta <- solve(H, -g)
    
    theta <- theta + delta
  }
  return(theta)
}

# Gauss-Newton Method (nls)
GN <- function(f, inits, data, minimum=TRUE, tol=1e-5){
  theta <- inits
  delta <- 10^3
  
  if (minimum == FALSE){
    f <- -f
  }
  
  while (min(abs(delta/theta)) > tol) {
    j <- jacobian(function(x, data) data$y - f(data$x, x), theta, data = data)
    r <- data$y - f(data$x, theta)
    
    g <- t(j) %*% r
    H <- t(j) %*% j
    
    delta <- solve(H, -g)
    
    theta <- theta + delta
  }
  return(theta)
}
