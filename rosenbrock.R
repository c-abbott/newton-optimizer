# rosenbrock function
rosenbrock <- function(theta, k){
  # split input vector to z and x so that the function works and the deriv function understands
  z <- theta[1]; x <- theta[2]
  f <- k * (z - x^2)^2 + (1 - x)^2
  ## create attributes for the f function and store the values for the gradient and the hessian matrix

  # attr(f, 'gradient') <- c(attr(Deriv.wrapper(theta,k), "gradient"))
  # attr(f, 'hessian') <- matrix(attr(Deriv.wrapper(theta,k), "hessian"),nrow = length(theta), ncol = length(theta))

  ## we can also use the following but i don't know which is more correct:

  attr(f, 'gradient') <- c(attr(grb(z,x,k), "gradient"))
  attr(f, 'hessian') <- matrix(attr(grb(z,x,k), "hessian"),nrow = length(theta), ncol = length(theta))

  ## we can use the above because the wrapper function does the same that happens in line 4

  ## activate for version where only gradient is returned and not the Hessian matrix

  # gradient calculated by hand

  # attr(f, 'gradient') <- c(k * (2 * (z - x^2)), -(2 * (1 - x) + k * (2 * (2 * x * (z - x^2)))))

  return(f)
}


# use of deriv function to create the derivatives and the Hessian matrix
# this is a function.
grb <- deriv(expression(k * (z - x^2)^2 + (1 - x)^2),c("z", "x"), ## dif wrt these
             function(z,x,k){}, hessian = TRUE)



# wrapper function for grb function so that we can use the z,x from theta.
# not needed probably.....
Deriv.wrapper <- function(theta, k) {
  z <- theta[1]; x <- theta[2]
  return(grb(z,x,k))
}


