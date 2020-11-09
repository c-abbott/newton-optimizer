
rosenbrock <- function(theta, k){
  z <- theta[1]; x <- theta[2]
  k * (z - x^2)^2 + (1 - x)^2
}

f = expression(k * (z - x^2)^2 + (1 - x)^2)

# deriv as an expression output
# grb <- deriv(expression(k*(z - x^2)^2 + (1 - x)^2),
#             c("z", "x"), ## dif wrt these
#             hessian = TRUE)

# this also works but grb is a function instead
# grb <- deriv(expression(k*(z - x^2)^2 + (1 - x)^2),
#              c("z", "x"), ## dif wrt these
#              function(z,x,k){}, hessian = TRUE)
# this does not need wrapper function for k, it is passed as argument

# wrapper function for deriv as a function
wrapper <- function(theta, k) {
  z <- theta[1]; x <- theta[2]
  grb <- deriv(f,
               c("z", "x"), ## dif wrt these
               function(z,x,k){}, hessian = TRUE)
  return(grb(z,x,k))
}




# # wrapper function for deriv as expression
# testf <- function(theta,k) {
#   z <- theta[1]; x <- theta[2] # need to code it to inherit from rosenbrock
#   mat <- eval(grb)
#   vec <- c(attr(mat,"gradient")) # create a vector with the gradient
#   # create a matrix with the Hessian
#   Hessian <- matrix(attr(mat, "hessian"), nrow = 2, ncol = 2) # rows and columns should be length or similar and dynamic
#   # R cannot output multiple objects so we create a list
#   output <- list("gradient" = vec,"Hessian" = Hessian )
#   return(output)
# }






# example of derivative as function
(fx <- deriv(y ~ b0 + b1 * 2^(-x/th), c("b0", "b1", "th"),function(b0, b1, th, x = 1:7){} ) )
fx(2, 3, 4)




newton <- function(theta, f, ..., tol = 1e-8, fscale = 1, maxit = 100, max.half = 20){

  while (n < maxit) {
    f0 <- rosenbrock(theta, k)
    g <- attr(wrapper(theta,k), "gradient")
    # create if not NULL statement
    Hs <- matrix(attr(wrapper(theta,k), "hessian"),nrow = length(theta), ncol = length(theta)) # rows and columns should be length or similar and dynamic
    if (abs(g) < tol) {
      break(0)
    }
    th1 <- th0 - (f0/g)
    if (max(abs(g)) < (abs(f0)+fscale)*tol){
      cat("converged")
      break(0)
    }
    th0 <- th1
    # also print iteration at this point
    print(th0)
    n <- n + 1
  }


  }

n = 1
fscale = 1
tol = 1e-8

# test while loop outside the newton function
while (n<maxit) {
  n <- n + 1
  f0 <- rosenbrock(th0,k)
  g <- attr(wrapper(th0,k), "gradient")
  # create if not NULL statement
  Hs <- matrix(attr(wrapper(th0,k), "hessian"),nrow = length(th0), ncol = length(th0)) # rows and columns should be length or similar and dynamic
  if (abs(g) < tol) {
    break(0)
  }
  th1 <- th0 - (f0/g) # objective function should be replaced
  if (max(abs(g)) < (abs(f0)+fscale)*tol){ #the condition has length > 1 and only the first element will be used
    cat("converged")
    break(0)
  }
  th0 <- th1
  # also print iteration at this point
  print(th0)
}




# lecture examples
## finite differencing
th0 <- c(-5,5)
fd.rb <- th0*0
rb.lo <- rosenbrock(th0,10)
eps <- 1e-8
for (i in 1:length(th0)){
  th1 <- th0;th1[i] <- th1[i]+eps
  rb.hi <- rosenbrock(th1,10)
  fd.rb[i] <- (rb.hi - rb.lo)/eps
}




grb0.1 <- function(theta,k) {
  z <- theta[1]; x <- theta[2]
  c(k * (2 * (z - x^2)), -(2 * (1 - x) + k * (2 * (2 * x * (z - x^2)))))
}


newton <- function(theta, f, ..., tol=1e-8, fscale=1, maxit=100, max.half=20) {
  # Start by assuming that the function passed by the user holds a gradient
  # and hessian attribute
  iter = 0
  while (iter < maxit) {
    # Getting gradient, hessian and inverse hessian
    gradient <- attr(f, 'gradient')
    H <- attr(f, 'hessian')
    Hi <- solve(H)

    # Evaluating objective function at initial guess
    f0 <- f(theta, ...)

    # Convergence check
    if (max(abs(gradient)) < (abs(f0)+fscale)*tol){
      cat("converged")
      return(list(f0, theta, iter, gradient, Hi))
    } else {
      # Checking hessian is positive definite
      is.pos.def = FALSE
      while (is.pos.def == FALSE) {
        tryCatch(
          expr <- {
            H.chol = chol(H)
            is.pos.def = TRUE
          },
          error <- function(e) {
            # Perturbing Hessian to be positive definite
            H <- H + diag(abs(max(H))*10*1e-9, nrow=nrow(H), ncol=ncol(H))
          }
        )
      }
      Delta <- backsolve(H.chol, forwardsolve(t(H.chol), -gradient))
      while (f(theta + Delta, ...) < f0) {
        Delta <- Delta / 2
      }
      theta = theta + Delta
    }
  }
}





