#Plate-Shaped

source('newton.R')


grb <- deriv(expression(k1 * (x1^2 + x2^2)+k2*x1*x2),c("x1", "x2"),
             function(x1,x2,k1,k2){}, hessian = TRUE)


matya1 <- function(theta, k1,k2){

  x1 <- theta[1]; x2 <- theta[2]
  f <- k1 * (x1^2 + x2^2)+k2*x1*x2

  attr(f, 'gradient') <- c(attr(grb(x1,x2,k1,k2), "gradient"))
  attr(f, 'hessian') <- matrix(attr(grb(x1,x2,k1,k2), "hessian"),nrow = length(theta), ncol = length(theta))
  return(f)
}


matya2 <- function(theta, k1,k2){

  x1 <- theta[1]; x2 <- theta[2]
  f <- k1 * (x1^2 + x2^2)+k2*x1*x2

  attr(f, 'gradient') <- c(attr(grb(x1,x2,k1,k2), "gradient"))
  return(f)
}


output1 <- newton(c(20, 20), matya1, k1=0.26,k2=0.48)
output2 <- newton(c(20, 20), matya2, k1=0.26,k2=0.48)


