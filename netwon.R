newton <- function(theta, f, ..., tol=1e-8, fscale=1, maxit=100, max.half=20) {
  # Evaluating objective function at initial guess
  f0 <- f(theta, ...)
  # Getting gradient, hessian and inverse hessian
  gradient <- attr(f0, 'gradient')
  H <- attr(f0, 'hessian')
  Hi <- solve(H)

  iter = 0
  while (iter < maxit) {
    # Convergence check
    if (max(abs(gradient)) < (abs(f0)+fscale)*tol){
      cat("Converged")
      return(list(f0, theta, iter, gradient, Hi))
    } else {
      # Checking hessian is positive definite
      is.pos.def = FALSE
      while (is.pos.def == FALSE) {
        H.chol <- tryCatch(
          {
            chol(H)
            is.pos.def = TRUE
            return(chol(H))
          },
          error <- function(cond) {
            # Perturbing Hessian to be positive definite
            H <- H + diag(abs(max(H))*1e-8*10^iter, nrow=nrow(H), ncol=ncol(H))
            return(H)
          }
        )
      }
      # Calculate forward step
      Delta <- backsolve(H.chol, forwardsolve(t(H.chol), -gradient))
      # Ensure steps are taken in right direction towards optimum
      while (f(theta + Delta, ...) < f0) {
        Delta <- Delta / 2
      }
      # Updating theta
      theta = theta + Delta
      iter <- iter + 1
    }
  }
}