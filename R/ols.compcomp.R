ols.compcomp <- function(y, x, xnew = NULL) {

  py <- dim(y)[2]   ;    px <- dim(x)[2]
  pyx <- py * px

  ols <- function(be) {
    be <- matrix(be, ncol = py)
    be <- be / rowsums(be)
    mu <- x %*% be
    sum( (y - mu)^2 )
  }

  runtime <- proc.time()
  mod <- optim( runif(pyx), ols, method = "L-BFGS-B", lower = rep(0, pyx),
                upper = rep(1, pyx), control = list(maxit = 10000) )
  mod <- optim( mod$par, ols, method = "L-BFGS-B", lower = rep(0, pyx),
                upper = rep(1, pyx), control = list(maxit = 10000) )
  runtime <- proc.time() - runtime

  be <- mod$par
  be <- matrix(be, ncol = py)
  be <- be / rowsums(be)


  if ( is.null( colnames(y) ) ) {
    colnames(be) <- paste("Y", 1:py, sep = "")
  } else colnames(be) <- colnames(y)
  if ( is.null( rownames(y) ) ) {
    rownames(be) <- paste("X", 1:px, sep = "")
  } else rownames(be) <- colnames(x)

  est <- NULL
  if ( !is.null(xnew) ) {
    est <- xnew %*% be
  }

  list( runtime = runtime, mse = mod$value / dim(y)[1], be = be, est = est )
}
