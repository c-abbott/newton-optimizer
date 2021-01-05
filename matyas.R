#Plate-Shaped

source('newton.R')


grb <- deriv(expression(k * (x1^2 + x2^2)+k*x1*x2),c("x1", "x2"),
             function(x1,x2,k){}, hessian = TRUE)


matya1 <- function(theta, k){
  x1 <- theta[1]; x2 <- theta[2]
  f <- k * (x1^2 + x2^2)+k*x1*x2
  attr(f, 'gradient') <- c(attr(grb(x1,x2,k), "gradient"))
  attr(f, 'hessian') <- matrix(attr(grb(x1,x2,k), "hessian"),nrow = length(theta), ncol = length(theta))
  return(f)
}


matya2 <- function(theta, k){
  x1 <- theta[1]; x2 <- theta[2]
  f <- k * (x1^2 + x2^2)+k*x1*x2

  attr(f, 'gradient') <- c(attr(grb(x1,x2,k), "gradient"))
  return(f)
}

output1 <- newton(c(20, 20), matya1, k=0.26)
output2 <- newton(c(20, 20), matya2, k=0.26)


