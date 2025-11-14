### First Draft of functions
library(numDeriv)

Newton <- function(start, f, data, epsilon){
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

Newton2D <- function(start, f, data, epsilon){
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

GaussNewton <- function(){
  
}