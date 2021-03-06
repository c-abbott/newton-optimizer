# Finite difference Hessian function (called when the f passed to newton() has
# no hessian attribute attached)
fd.Hessian <- function(theta, f, k) {
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
  if (any(is.null(H))) {
    # Calling fd.Hessian defined above
    H <- fd.Hessian(theta, f,...)
  }

  # Checking quantities needed for optimization are finite
  if (!is.finite(f0) || any(!is.finite(gradient))|| any(!is.finite(H))) {
    stop("The objective or its derivatives are not finite at
         your initial estimate of theta.")
  }

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
      # Inverse Hessian calculation
      Hi <- solve(H)
      return(list(f.min=f0[1], theta=theta, iter=iter, g=gradient, Hi=Hi))
    } else {
      catch.iter <- 0
      # Cholesky factor calculation (perturb H until this calculation is possible)
      while(inherits(try(H.chol <- chol(H), silent = TRUE), "try-error") == TRUE) {
        # Perturbing hessian to be positive definite (runs of if chol(H) fails)
        H <- H + diag(abs(max(H))*10^catch.iter*1e-8, nrow=nrow(H), ncol=ncol(H))
        catch.iter <- catch.iter + 1
      }
      # Calculate forward step towards optimum (minimum)
      Delta <- backsolve(H.chol, forwardsolve(t(H.chol), -gradient))

      # Ensuring steps are taken in right direction towards optimum (minimum)
      # a.k.a preventing optimizer from blowing up
      half.iter = 0
      while ((f(theta + Delta, ...)[1] > f0[1]) |
             !is.finite(f(theta + Delta, ...)) |
             any(!is.finite(attr(f(theta + Delta, ...), "gradient"))) |
             any(!is.finite(attr(f(theta + Delta, ...), "hessian")))) {
        if (half.iter < max.half) {
          Delta <- Delta / 2
          half.iter <- half.iter + 1
        } else {
          stop(paste("The update step failed to reduce the objective
                  after ", as.character(max.half), " halvings"))
        }
      }
      # Updating theta
      theta <- theta + Delta

      # Updating function values and attributes at new theta
      f0 <- f(theta, ...)
      gradient <- attr(f0, 'gradient')

      # Hessian calculation by finite differencing if f has no hessian attribute
      if (any(is.null(attr(f0, 'hessian')))) {
        # Calling fd.Hessian defined above
        H <- fd.Hessian(theta, f,...)
      } else {
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