#Steep Ridges/Drops
source('newton.R')

grb <- deriv(expression(cos(x)*cos(y)*exp(-(x-k)^2-(y-k)^2)),c("x", "y"),
             function(x,y,k){}, hessian = TRUE)

easom1 <- function(theta,k){

  x <- theta[1];y <- theta[2]

  f <- cos(x)*cos(y)*exp(-(x-k)^2-(y-k)^2)

  attr(f, 'gradient') <- c(attr(grb(x,y,k), "gradient"))
  attr(f, 'hessian') <- matrix(attr(grb(x,y,k), "hessian"),nrow = length(theta), ncol = length(theta))
  return(f)
}

easom2 <- function(theta,k){

  x <- theta[1]; y <- theta[2]

  f <- cos(x)*cos(y)*exp(-(x-k)^2-(y-k)^2)

  attr(f, 'gradient') <- c(attr(grb(x,y,k), "gradient"))

  return(f)
}

# Calls to newton optimizer
output1 <- newton(c(20, 10), easom1, k = 3)
output2 <- newton(c(20, 10), easom2, k = 3)

