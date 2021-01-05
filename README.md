[![Check Required Files](https://github.com/statprog-s1-2020/hw03_tut01_team08/workflows/Check%20Required%20Files/badge.svg)](https://github.com/statprog-s1-2020/hw03_tut01_team08/actions?query=workflow:%22Check%20Required%20Files%22) [![Check Scripts Run](https://github.com/statprog-s1-2020/hw03_tut01_team08/workflows/Check%20Scripts%20Run/badge.svg)](https://github.com/statprog-s1-2020/hw03_tut01_team08/actions?query=workflow:%22Check%20Scripts%20Run%22)


Statistical Programming 
---------

* Christos Konstantinou 
* Callum Abbott 
* Zhang Yuying 
* Guanxing Huang 


## Writing a Newton optimizer

### Task:

Write an R function, `newton`, implementing Newton’s method for minimization of functions, and provide example code using it to optimize Rosenbrock’s function.

### Specification: 

Your `newton` optimization function should operate broadly in the same way as `nlm`. Note that the purpose is to have an independent implementation: you must code the optimization yourself, not simply call optimization code written by someone else. 

The function signature should be as follows,
```r
newton(theta, f, ..., tol = 1e-8, fscale = 1, maxit = 100, max.half = 20)
```
with the arguments defined as follows:

* `theta` is a vector of initial values for the optimization parameters.
* `f` is the objective function to minimize. Its first argument is the vector of optimization parameters. Remaining arguments will be passed from newton using `...`. The scalar value returned by `f` will have two attributes, a `"gradient"` vector and, optionally, a `"hessian"` matrix.
* `...` any arguments of `f` after the first (the parameter vector) are passed using this (see hints below). 
* `tol` the convergence tolerance.
* `fscale` a rough estimate of the magnitude of `f` at the optimum - used in convergence testing.
* `maxit` the maximum number of Newton iterations to try before giving up.
* `max.half` the maximum number of times a step should be halved before concluding that the step has failed to improve the objective.


Your `newton` function should return a list containing:
* `f` the value of the objective function at the minimum.
* `theta` the value of the parameters at the minimum.
* `iter` the number of iterations taken to reach the minimum.
* `g` the gradient vector at the minimum (so the user can judge closeness to  numerical zero).
* `Hi` the inverse of the Hessian matrix at the minimum (useful if the objective is a negative log likelihood).

The function should issue errors or warnings (using `stop` or `warning` as appropriate) in at least the following cases. 

1. If the objective or derivatives are not finite at the initial `theta`. 
2. If the step fails to reduce the objective despite trying `max.half` step halvings. 
3. If `maxit` is reached without convergence.
4. If the Hessian is not positive definite at convergence.


### Other considerations:

To judge whether the gradient vector is close enough to zero, you will need to consider the magnitude of the objective (you can't expect gradients to to be down at 10<sup>-10</sup> if the objective is of order 10<sup>10</sup>, for example). So the gradients are judged to be zero when they are smaller than `tol` multiplied by the objective. But then there is a problem it the objective is zero at the minimum - we can never succeed in making the magnitude of the gradient less than the magnitude of the objective. So `fscale` is provided. Then if `f0` is the current value of the objective and `g` the current gradient,
    
    ```
    max(abs(g)) < (abs(f0)+fscale)*tol
    ``` 
    
    is a suitable condition for convergence.

If no Hessian matrix is supplied, your code should generate one by finite differencing the gradient vector. Such an approximate Hessian will be asymmetric: `H <- 0.5 * (t(H) + H)` fixes that.

