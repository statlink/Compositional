tflr.irls <- function(y, x, xnew = NULL, tol = 1e-08, maxit = 100) {

  runtime <- proc.time()
  py <- dim(y)[2]   ;    px <- dim(x)[2]
  pyx <- py * px    ;    n <- dim(y)[1]

  dvec <- as.vector( crossprod(x, y) )
  xx <- crossprod(x)

  XX <- matrix(0, pyx, pyx)
  ind <- matrix( 1:pyx, ncol = px, byrow = TRUE )
  for ( i in 1:py )  XX[ ind[i, ], ind[i, ] ] <- xx
  A <- matrix(0, pyx, pyx)
  for ( i in 1:px )  A[i, ind[, i]] <- 1
  A <- t( rbind( A, diag(pyx) ) )
  A <- A[, -c( (px + 1): pyx) ]
  bvec <- c( rep(1, px), rep(0, pyx) )
  f <- try( quadprog::solve.QP( Dmat = XX, dvec = dvec, Amat = A, bvec = bvec,
                                meq = px ), silent = TRUE )
  if ( identical(class(f), "try-error") ) {
    f <- quadprog::solve.QP( Dmat = Matrix::nearPD(XX)$mat, dvec = dvec, Amat = A, bvec = bvec, meq = px )
  }
  be1 <- matrix( abs(f$solution), ncol = py)

  pi_hat <- x %*% be1
  #dev1 <- sum( y * log(y / pi_hat), na.rm = TRUE )
  a <- y * log(pi_hat)
  a[is.infinite(a)] <- NA
  dev1 <-  - sum(a, na.rm = TRUE)

  # Working weights and response
  w <- 1 / (pi_hat * (1 - pi_hat) )
  # z <- y
  # Weighted least squares update: beta = (X'WX)^(-1) X'Wz
  #wx <- w * x
  for ( i in 1:py )  XX[ ind[i, ], ind[i, ] ] <- crossprod(x, w[, i] * x)
  dvec <- as.vector( crossprod( x, w * y ) )
  f <- try( quadprog::solve.QP( Dmat = XX, dvec = dvec, Amat = A, bvec = bvec,
                                meq = px ), silent = TRUE )
  if ( identical(class(f), "try-error") ) {
    f <- quadprog::solve.QP( Dmat = Matrix::nearPD(XX)$mat, dvec = dvec, Amat = A, bvec = bvec, meq = px )
  }
  be2 <- matrix( f$solution, ncol = py)
  pi_hat <- x %*% be2
  #dev2 <- sum( y * log(y / pi_hat) )
  a <- y * log(pi_hat)
  a[is.infinite(a)] <- NA
  dev2 <-  - sum(a, na.rm = TRUE)
  i <- 2

  # IRLS iteration
  while ( dev1 - dev2 > tol  &  i < maxit ) {
    i <- i + 1
    dev1 <- dev2
    be1 <- be2
    w <- 1 / ( pi_hat * (1 - pi_hat) )
    for ( i in 1:py )  XX[ ind[i, ], ind[i, ] ] <- crossprod(x, w[, i] * x)
    dvec <- as.vector( crossprod( x, w * y ) )
    f <- try( quadprog::solve.QP( Dmat = XX, dvec = dvec, Amat = A, bvec = bvec,
                                  meq = px ), silent = TRUE )
    if ( identical(class(f), "try-error") ) {
      f <- quadprog::solve.QP( Dmat = Matrix::nearPD(XX)$mat, dvec = dvec, Amat = A, bvec = bvec, meq = px )
    }
    be2 <- matrix( abs(f$solution), ncol = py)
    pi_hat <- x %*% be2
    #dev2 <- sum( y * log(y / pi_hat) )
    a <- y * log(pi_hat)
    a[is.infinite(a)] <- NA
    dev2 <-  - sum(a, na.rm = TRUE)
  }
  runtime <- proc.time() - runtime

  if ( dev1 > dev2 ) {
    be <- matrix(be2, ncol = py)
  } else  be <- matrix(be1, ncol = py)

  if ( is.null( colnames(y) ) )  {
    colnames(be) <- paste("Y", 1:py, sep = "")
  } else colnames(be) <- colnames(y)

  if ( is.null( colnames(x) ) )  {
    rownames(be) <- paste("X", 1:px, sep = "")
  } else rownames(be) <- colnames(x)
  est <- NULL
  if ( !is.null(xnew) )  est <- xnew %*% be
  kl <- sum( y * log(y), na.rm = TRUE ) + dev2

  list(runtime = runtime, iters = i, kl = kl, be = be, est = est)
}
