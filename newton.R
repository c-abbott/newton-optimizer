# Finite difference Hessian function (called when the f passed to newton() has
# no hessian attribute attached)
fd.Hessian <- function(theta,f,k) {
  eps <- 1e-8
  gradient <- attr(f(theta, k), 'gradient')
  # Initialize Hessian matrix
  fd.f <- matrix(0L, nrow = length(theta), ncol = length(theta))
  for (i in 1:length(theta)){
    th1 <- theta; th1[i] <- th1[i]+eps
    f.hi <- attr(f(th1,k), 'gradient')
    fd.f[,i] <- (f.hi - gradient)/eps
  }
  # Make created Hessian matrix symmetric
  fd.f <- 0.5 * (t(fd.f) + fd.f)
  return(fd.f)
}

newton <- function(theta, f, ..., tol=1e-8, fscale=1, maxit=100, max.half=20) {
  # Evaluating objective function and attributes at initial theta
  f0 <- f(theta, ...)
  gradient <- attr(f0, 'gradient')
  H <- attr(f0, 'hessian')

  # Hessian calculation by finite differencing if f has no hessian attribute
  if (is.null(H)) {
    # Calling fd.Hessian defined above
    attr(f0, 'hessian') <- fd.Hessian(theta, f,...)
    H <- attr(f0, 'hessian')
  }

  # Checking quantities needed for optimization are finite
  if (!is.finite(f0) || !is.finite(gradient) || !is.finite(H)) {
    stop("The objective or its derivatives are not finite at
         your initial estimate of theta.")
  }
  # Inverse Hessian calculation
  Hi <- solve(H)

  # ------------------- Optimization Begins ------------------- #
  iter = 0
  while (iter < maxit) {
    # Convergence check
    if (max(abs(gradient)) < (abs(f0)+fscale)*tol){
      cat("Converged \n")
      # Checking hessian is positive definite at convergence
      tryCatch({
        chol(H)
      },
      error = function(cond){
        warning("The Hessian is not positive definite at convergence")
      })
      return(list(f0, theta, iter, gradient, Hi))
    } else {
      # Checking hessian is positive definite by observing whether its
      # Cholesky factor exists
      is.pos.def = FALSE
      while (is.pos.def == FALSE) {
        H <- tryCatch({
          chol(H)
          is.pos.def = TRUE
          H
        },
        error = function(cond){
          # Perturbing Hessian to be positive definite if chol(H) does not exist
          H <- H + diag(abs(max(H))*1e-8*10*iter, nrow=nrow(H), ncol=ncol(H))
        })
      }
      # Calculate forward step towards optimum (minimum)
      H.chol <- chol(H)
      Delta <- backsolve(H.chol, forwardsolve(t(H.chol), -gradient))

      # Ensuring steps are taken in right direction towards optimum (minimum)
      # a.k.a preventing optimizer from blowing up
      half.iter = 0
      while (f(theta + Delta, ...) > f0) {
        if (half.iter < max.half) {
          Delta <- Delta / 2
          half.iter <- half.iter + 1
        } else {
          warning(paste("The update step failed to reduce the objective
                  after ", as.character(max.half), " halvings"))
        }
      }
      # Updating theta
      theta <- theta + Delta

      # Updating function values and attributes at new theta
      f0 <- f(theta, ...)
      gradient <- attr(f0, 'gradient')
      H <- attr(f0, 'hessian')
      # Hessian calculation by finite differencing if f has no hessian attribute
      if (is.null(H)) {
        # Calling fd.Hessian defined above
        attr(f0, 'hessian') <- fd.Hessian(theta, f,...)
        H <- attr(f0, 'hessian')
      }
      # Iteration update
      iter <- iter + 1
    }
  }
  # Checking whether convergence occurs when iter == maxit
  if (max(abs(gradient)) < (abs(f0)+fscale)*tol){
    cat("Converged")
    return(list(f0, theta, iter, gradient, Hi))
  } else {
    warning(paste("Newton optimizer failed to converage after
                  maxit = ", as.character(maxit), " iterations"))
  }
}