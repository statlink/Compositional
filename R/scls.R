scls <- function(y, x, xnew = NULL, nbcores = 4) {

  if ( inherits(y, "sparseMatrix")  |  inherits(x, "sparseMatrix") ) {
    res <- .scls_osqp(y, x, xnew)

  } else {
    py <- dim(y)[2]   ;    px <- dim(x)[2]
    pyx <- py * px    ;    n <- dim(y)[1]

    if ( identical( class(x)[1], "FBM") ) {
      xx <- bigstatsr::big_crossprodSelf(x)
      xx <- xx[1:px, 1:px]
      a1 <- bigstatsr::FBM(px, py)
      dvec <- bigstatsr::big_apply(y, function(X, ind, x, res) {
        res[, ind] <- bigstatsr::big_cprodMat(x, X[, ind, drop = FALSE])
      }, a.combine  = "c", block.size = 500, ncores = nbcores, x = x, res = a1)
      dvec <- dvec[1:pyx]
    } else {
      dvec <- as.vector( crossprod(x, y) )
      xx <- crossprod(x)
    }

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

    be <- matrix( abs(f$solution), ncol = py)
    mse <- ( sum(y^2) + 2 * f$value ) / n

    if ( is.null( colnames(y) ) ) {
      colnames(be) <- paste("Y", 1:py, sep = "")
    } else colnames(be) <- colnames(y)
    if ( is.null( colnames(x) ) ) {
      rownames(be) <- paste("X", 1:px, sep = "")
    } else rownames(be) <- colnames(x)

    est <- NULL
    if ( !is.null(xnew) ) {
      est <- xnew %*% be
    }

    res <- list( mse = mse, be = be, est = est )
  }
  res
}




.scls_osqp <- function(y, x, xnew = NULL) {

  py <- ncol(y)    ;    px <- ncol(x)
  pyx <- py * px   ;    n <- nrow(y)

  # Convert inputs to sparse if not already
  if ( !inherits(y, "sparseMatrix") )  y <- Matrix::Matrix(y, sparse = TRUE)
  if ( !inherits(x, "sparseMatrix") )  x <- Matrix::Matrix(x, sparse = TRUE)

  xx <- Matrix::crossprod(x)  # sparse
  xy <- Matrix::crossprod(x, y)  # sparse

  # Construct sparse P via Kronecker product
  P <- 2 * Matrix::kronecker(Matrix::Diagonal(py), xx)
  # q vector (usually dense, but keep efficient)
  q <- -2 * as.vector( as.matrix(xy) )
  # Equality: sum-to-one for each predictor
  A_eq <- Matrix::kronecker(Matrix::Diagonal(px), Matrix::Matrix(rep(1, py), nrow = 1, sparse = TRUE))
  # Inequality: non-negativity
  A_ineq <- Matrix::Diagonal(pyx)
  # Stack constraints (both sparse)
  A <- rbind(A_eq, A_ineq)
  l <- c(rep(1, px), rep(0, pyx))
  u <- c(rep(1, px), rep(Inf, pyx))
  # Solve with OSQP
  settings <- osqp::osqpSettings(verbose = FALSE, eps_abs = 1e-8, eps_rel = 1e-8)
  model <- osqp::osqp(P = P, q = q, A = A, l = l, u = u, pars = settings)
  res <- model$Solve()
  be <- matrix(res$x, ncol = py)
  #fitted <- x %*% be
  #mse <- sum( (as.matrix(y) - as.matrix(fitted))^2 ) / n
  mse <- ( model$Solve()$info$obj_val + sum(y^2) ) / n

  if (is.null(colnames(y))) {
    colnames(be) <- paste0("Y", 1:py)
  } else colnames(be) <- colnames(y)

  if (is.null(colnames(x))) {
    rownames(be) <- paste0("X", 1:px)
  } else rownames(be) <- colnames(x)

  est <- NULL
  if ( !is.null(xnew) ) {
    if (!inherits(xnew, "sparseMatrix")) {
      xnew <- Matrix::Matrix(xnew, sparse = TRUE)
    }
    est <- as.matrix(xnew %*% be)
  }

  list(mse = mse, be = be, est = est)
}







# ols.compcomp_newer <- function(y, x, xnew = NULL) {

  # py <- dim(y)[2]   ;    px <- dim(x)[2]
  # pyx <- py * px

  # xxinv <- solve( crossprod(x) )
  # m <- numeric( px )
  # dvec <-  - 2 * as.vector( crossprod(x, y) )
  # yy <- sum( diag( crossprod(y) ))

  # ols <- function(be, yy, dvec, m, xxinv) {
    # be <- matrix(be, ncol = py)
    # be <- be / rowSums(be)
    # be1 <- as.vector(be)
    # yy + sum(dvec * be1) + sum( Rfast::mahala(t(be), m, xxinv) )
  # }

  # runtime <- proc.time()
  # mod <- optim( runif(pyx), ols, yy = yy, dvec = dvec, m = m, xxinv = xxinv,
                # method = "L-BFGS-B", lower = rep(0, pyx),
                # upper = rep(1, pyx), control = list(maxit = 10000) )
  # mod <- optim( mod$par, ols, yy = yy, dvec = dvec, m = m, xxinv = xxinv,
                # method = "L-BFGS-B", lower = rep(0, pyx),
                # upper = rep(1, pyx), control = list(maxit = 10000) )
  # runtime <- proc.time() - runtime

  # be <- mod$par
  # be <- matrix(be, ncol = py)
  # be <- be / rowsums(be)


  # if ( is.null( colnames(y) ) ) {
    # colnames(be) <- paste("Y", 1:py, sep = "")
  # } else colnames(be) <- colnames(y)
  # if ( is.null( rownames(y) ) ) {
    # rownames(be) <- paste("X", 1:px, sep = "")
  # } else rownames(be) <- colnames(x)

  # est <- NULL
  # if ( !is.null(xnew) ) {
    # est <- xnew %*% be
  # }

  # list( runtime = runtime, mse = mod$value / dim(y)[1], be = be, est = est )
# }




# ols.compcomp_old <- function(y, x, xnew = NULL) {

  # py <- dim(y)[2]   ;    px <- dim(x)[2]
  # pyx <- py * px

  # ols <- function(be) {
    # be <- matrix(be, ncol = py)
    # be <- be / rowsums(be)
    # mu <- x %*% be
    # sum( (y - mu)^2 )
  # }

  # runtime <- proc.time()
  # mod <- optim( runif(pyx), ols, method = "L-BFGS-B", lower = rep(0, pyx),
                # upper = rep(1, pyx), control = list(maxit = 10000) )
  # mod <- optim( mod$par, ols, method = "L-BFGS-B", lower = rep(0, pyx),
                # upper = rep(1, pyx), control = list(maxit = 10000) )
  # runtime <- proc.time() - runtime

  # be <- mod$par
  # be <- matrix(be, ncol = py)
  # be <- be / rowsums(be)


  # if ( is.null( colnames(y) ) ) {
    # colnames(be) <- paste("Y", 1:py, sep = "")
  # } else colnames(be) <- colnames(y)
  # if ( is.null( rownames(y) ) ) {
    # rownames(be) <- paste("X", 1:px, sep = "")
  # } else rownames(be) <- colnames(x)

  # est <- NULL
  # if ( !is.null(xnew) ) {
    # est <- xnew %*% be
  # }

  # list( runtime = runtime, mse = mod$value / dim(y)[1], be = be, est = est )
# }
