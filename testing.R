### First Draft of functions
library(numDeriv)

# Univariate Newton Method
UVN <- function(f, inits, data, minimum=TRUE, tol=1e-5, maxit=10000){
  theta <- inits
  delta <- 10^3
  iter <- 0
  conv <- 0
  
  if (minimum == FALSE){
    f <- -f
  }
  
  while (abs(delta/theta) > tol && iter <= maxit){ #TRY/EXCEPT => 1
    g <- grad(f, theta, data = data)
    h <- hessian(f, theta, data = data)[1,1]
    
    delta <- -g/h
    theta <- theta + delta
    
    iter <- iter + 1
  }
  
  if (iter <= maxit){
    conv <- 2
  }
  
  return(list(estimate=theta, feval=f(theta, data$x), tolerance=abs(delta/theta),
              conv=conv, niter=iter))
}

# Multivariate Newton Method
MVN <- function(f, inits, data, minimum=TRUE, tol=1e-5, maxit=10000){
  theta <- inits
  delta <- 10^3
  iter <- 0
  conv <- 0
  
  if (minimum == FALSE){
    f <- -f
  }
  
  # For each step, calculate gradient, hessian, solve linear system for delta
  while (min(abs(delta/theta)) > tol && iter <= maxit){ # Take minimum abs. change in each 
    g <- grad(f, theta, data = data)       # par, not total step size
    H <- hessian(f, theta, data = data)
    
    delta <- solve(H, -g)
    theta <- theta + delta
    
    iter <- iter + 1
  }
  
  if (iter <= maxit){
    conv <- 2
  }
  
  return(list(estimate=theta, feval=f(theta, data$x), tolerance=abs(delta/theta),
              conv=conv, niter=iter))
}

# Gauss-Newton Method (nls)
GN <- function(f, inits, data, minimum=TRUE, tol=1e-5, maxit=10000){
  theta <- inits
  delta <- 10^3
  iter <- 0
  conv <- 0
  
  if (minimum == FALSE){
    f <- -f
  }
  
  while (min(abs(delta/theta)) > tol && iter <= maxit) {
    j <- jacobian(f, theta, data = data) #function(x, data) data$y - f(data$x, x)
    r <- data$y - f(data$x, theta)
    
    g <- t(j) %*% r
    H <- t(j) %*% j
    
    delta <- solve(H, -g)
    theta <- theta + delta
    
    iter <- iter + 1
  }
  
  if (iter <= maxit){
    conv <- 2
  }
  
  return(list(estimate=theta, feval=f(theta, data$x), tolerance=abs(delta/theta),
              conv=conv, niter=iter))
}
