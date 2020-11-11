source('newton.R')

grb <- deriv(expression(k*x1^2 -1.05*x1^4 + x1^6/6 + x1*x2 + x2^2,c("x1", "x2"),
             function(x1,x2,k){}, hessian = TRUE)


#Valley-Shaped
came <- function(theta,k){


x1 <- theta[1]; x2 <- theta[2]



f <- k*x1^2 -1.05*x1^4 + x1^6/6 + x1*x2 + x2^2



attr(f, 'gradient') <- c(attr(grb(x1,x2,k), "gradient"))
attr(f, 'hessian') <- matrix(attr(grb(x1,x2,k), "hessian"),nrow = length(theta), ncol = length(theta))
return(f)
}



came2 <- function(theta,k){



x1 <- theta[1]; x2 <- theta[2]



f <- k*x1^2 -1.05*x1^4 + x1^6/6 + x1*x2 + x2^2



attr(f, 'gradient') <- c(attr(grb(x1,x2,k), "gradient"))



return(f)
}



# Calls to newton optimizer
output1 <- newton(c(-1, 1), came, k = 2)
output2 <- newton(c(-1, 1), came2, k = 2)
