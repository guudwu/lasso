# Weighted LASSO solver via FISTA algorithm

fista_lasso <- function (
  y , x
  , beta = numeric(ncol(x))
  , w = 1
  , maxiter = 100
  , mu = svd(x)$d[1]**2
  , tol = 1e-5
)

# Solve weighted LASSO problem
#   \| y - x beta \|^2 + w^T abs(beta)
# by FISTA algorithm.

# INPUT:
# y: Dependent variable.
# x: Each column is an explanatory variable.
# beta: Initial guess of regression coefficients.
# w: Weight vector.
#   Might be a scalar, which indicates same weight on each coefficient.
# maxiter: Force to terminate before convergence
#   when number of iteration reaches this number.
# mu: Coefficient of quadratic term in approximation model.
#   Generally this controls the step-length of each iteration.
#   The larger mu is, the smaller the step-length.
#   Default value is the smallest value to guarantee convergence.
#   In practice, user may choose to further reduce it.
# tol: Precision tolerance.
#   Algorithm terminates when difference of objective function value
#   between adjacent steps is smaller than this number.

# OUTPUT:
# beta: Regression coefficients when algorithm terminates.
# fun: Objective function value when algorithm terminates.
#   First element is residual sum-of-squares.
#   Second element is LASSO term.
# iter: Number of iterations when algorithm terminates.

{
# Sanity check
  x <- as.matrix(x)

  y <- as.numeric(y)
  if ( length(y) != nrow(x) )
    stop('Dimension mismatch between argument "x","y".')

  beta <- as.numeric(beta)
  if ( length(beta) != ncol(x) )
    stop('Dimension mismatch between argument "x","beta".')

  w <- as.numeric(w)
  if ( length(w) == 1 )
    w <- rep(w,length(beta))
  if ( length(w) != length(beta) )
    stop('Dimension mismatch between argument "w","beta".')

  maxiter <- as.integer(maxiter)
  if ( length(maxiter) != 1 )
    stop('Argument "maxiter" must be an integer.')

  mu <- as.numeric(mu)
  if ( length(mu) != 1 )
    stop('Argument "mu" must be a scalar.')

  tol <- as.numeric(tol)
  if ( length(tol) != 1 )
    stop('Argument "tol" must be a scalar.')

# Objective function
  .objective <- function ( beta )
  {
    temp <- x %*% beta - y
    return( c ( sum(temp**2) , sum(abs(w*beta)) ) )
  }

  xtx <- t(x) %*% x
  xty <- t(x) %*% y

  objective_old <- .objective(beta)

  for ( count in 1 : maxiter )
  {
    grad_half <- xtx %*% beta - xty

    temp <- beta - grad_half/mu

# Shrinkage
    sign_temp <- sign(temp)
    temp <- abs(temp) - w/mu/2
    temp[which(temp<0)] <- 0
    beta <- sign_temp * temp

# Termination
    objective_new <- .objective(beta)
    if ( abs(sum(objective_new)-sum(objective_old)) < tol )
      break
    objective_old <- objective_new
  }

  return ( list (
    beta = beta
    , fun = objective_new
    , iter = count
  ) )
}
