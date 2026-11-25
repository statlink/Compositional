################################
#### Jensen-Shannon divergence based regression for compositional data
#### Tsagris Michail 5/2015
#### Newton-Raphson (analytic gradient and Hessian) version
#### mtsagris@yahoo.gr
#### References: Michail Tsagris (2015)
#### A novel, divergence based, regression for compositional data
#### Proceedings of the 28th Panhellenic Statistics Conference
################################
js.compreg <- function(y, x, con = TRUE, B = 1, ncores = 1, xnew = NULL,
                       tol = 1e-07, maxiters = 100) {

  ## Newton-Raphson fit; y: n x D, x: n x p design matrix, be: p x d starting values
  jsnr <- function(y, x, be, tol = 1e-07, maxiters = 100) {
    n <- dim(y)[1]
    p <- dim(x)[2]
    d <- dim(y)[2] - 1
    ylogy <- sum( y[y > 0] * log(y[y > 0]) )

    fit <- function(be) {
      eta <- cbind(0, x %*% be)
      eta <- eta - eta[ cbind(1:n, max.col(eta, ties.method = "first")) ]
      mu <- exp(eta)
      mu / rowSums(mu)
    }
    obj <- function(mu) {
      M <- 0.5 * (mu + y)
      ylogy - sum( y * log(M) ) + sum( mu * log(mu / M) )
    }
    hess <- function(m, A, gbar, abar) {
      H <- matrix(0, p * d, p * d)
      for (k in 1:d) {
        ik <- (k - 1) * p + 1:p
        for (l in k:d) {
          il <- (l - 1) * p + 1:p
          w <- - m[, k] * m[, l] * ( A[, k] + A[, l] - gbar - abar )
          if (k == l)  w <- w + m[, k] * ( A[, k] - gbar )
          blk <- crossprod(x, x * w)
          H[ik, il] <- blk
          if (k != l)  H[il, ik] <- t(blk)
        }
      }
      H
    }

    mu <- fit(be)
    f1 <- obj(mu)
    f0 <- f1 + 1
    i <- 0
    while ( abs(f0 - f1) > tol  &  i < maxiters ) {
      i <- i + 1
      f0 <- f1
      G <- log( 2 * mu / (mu + y) )
      C <- y / (mu + y)                ## = mu * h
      gbar <- rowSums(mu * G)
      cbar <- rowSums(mu * C)
      m <- mu[, -1, drop = FALSE]
      grad <- as.vector( crossprod( x, m * ( G[, -1, drop = FALSE] - gbar ) ) )
      H <- hess(m, G[, -1, drop = FALSE] + C[, -1, drop = FALSE], gbar, gbar + cbar)
      R <- tryCatch( chol(H), error = function(e) NULL )
      if ( is.null(R) ) {              ## Gauss-Newton fallback
        H <- hess(m, C[, -1, drop = FALSE], 0, cbar)
        R <- chol( H + diag(1e-10, p * d) )
      }
      step <- matrix( backsolve( R, forwardsolve( t(R), grad ) ), ncol = d )
      ## step-halving
      a <- 1
      repeat {
        benew <- be - a * step
        munew <- fit(benew)
        fnew <- obj(munew)
        if ( fnew <= f1  ||  a < 1e-10 )  break
        a <- a / 2
      }
      be <- benew
      mu <- munew
      f1 <- fnew
    }
    be
  }

  n <- dim(y)[1]  ## sample size
  x <- model.matrix(y ~ ., data.frame(x) )
  if ( !con )  x <- x[, -1, drop = FALSE]
  p <- dim(x)[2]
  d <- dim(y)[2] - 1  ## dimensionality of the simplex
  namx <- colnames(x)
  namy <- colnames(y)
  if ( is.null( namy ) )  {
    namy <- paste("Y", 2:(d + 1), sep = "")
  } else namy <- namy[-1]

  runtime <- proc.time()
  ## starting values from the KL regression
  if ( con ) {
    ini <- Compositional::kl.compreg(y, x[, -1, drop = FALSE], con = TRUE)$be
  } else  ini <- Compositional::kl.compreg(y, x, con = FALSE)$be
  ini <- matrix(ini, ncol = d)
  be <- jsnr(y, x, ini, tol = tol, maxiters = maxiters)
  covb <- NULL
  runtime <- proc.time() - runtime

  if (B > 1) {
    runtime <- proc.time()
    ## betas are stored row-wise, i.e. as.vector( t(be) ), to match the names below
    if (ncores <= 1) {
      betaboot <- matrix( nrow = B, ncol = p * d )
      for (i in 1:B) {
        ida <- rangen::Sample.int(n, n, replace = TRUE)
        bb <- jsnr(y[ida, ], x[ida, , drop = FALSE], be, tol = tol, maxiters = maxiters)
        betaboot[i, ] <- as.vector( t(bb) )
      }  ##  end for (i in 1:B) {
    } else {
      cl <- parallel::makeCluster(ncores)
      on.exit( parallel::stopCluster(cl), add = TRUE )
      parallel::clusterEvalQ(cl, library(rangen))
      parallel::clusterExport(cl, varlist = c("jsnr", "y", "x", "n", "be", "tol", "maxiters"),
                              envir = environment())
      betaboot <- t( parallel::parSapply(cl, 1:B, function(i) {
        ida <- rangen::Sample.int(n, n, replace = TRUE)
        bb <- jsnr(y[ida, ], x[ida, , drop = FALSE], be, tol = tol, maxiters = maxiters)
        as.vector( t(bb) )
      }) )
    }  ##  end if (ncores <= 1) {
    covb <- cov(betaboot)
    runtime <- proc.time() - runtime
    nam <- NULL
    for (i in 1:p)  nam <- c(nam, paste(namy, ":", namx[i], sep = "") )
    colnames(covb) <- rownames(covb) <- nam
  }  ##  end if (B > 1) {

  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )
    if ( !con )  xnew <- xnew[, -1, drop = FALSE]
    mu <- cbind( 1, exp(xnew %*% be) )
    est <- mu / Rfast::rowsums(mu)
    colnames(est) <- colnames(y)
  }  else  est <- NULL
  colnames(be) <- namy
  rownames(be) <- namx
  list(runtime = runtime, be = be, covbe = covb, est = est)
}









# js.compreg <- function(y, x, con = TRUE, B = 1, ncores = 1, xnew = NULL) {
#   ## y is dependent variable, the compositional data
#   ## x is the independent variable(s)
#   ## B is the number of bootstrap samples used to obtain
#   ## standard errors for the betas
#   ## if B==1 no bootstrap is performed and no standard errors are reported
#   ## if ncores=1, then 1 processor is used, otherwise
#   ## more are used (parallel computing)
#   jsreg <- function(para, y, x, d) {
#     be <- matrix(para, byrow = TRUE, ncol = d)
#     mu1 <- cbind( 1, exp(x %*% be) )
#     mu <- mu1 / rowSums(mu1)
#     M <- ( mu + y ) / 2
#     a <-  - y * log1p(mu / y)
#     a[y == 0] <- 0
#     sum( a + mu * log(mu / M) )
#   }
#
#   n <- dim(y)[1]  ## sample size
#   x <- model.matrix(y ~ ., data.frame(x) )
#   if ( !con )  x <- x[, -1, drop = FALSE]
#   p <- dim(x)[2]
#   d <- dim(y)[2] - 1  ## dimensionality of the simplex
#   namx <- colnames(x)
#   namy <- colnames(y)
#   if ( is.null( namy ) )  {
#     namy <- paste("Y", 2:(d + 1), sep = "")
#   } else namy <- namy[-1]
#   ## the next lines minimize the kl.compreg function and obtain the estimated betas
#   runtime <- proc.time()
#   ini <- as.vector( t( Compositional::kl.compreg(y, x[, -1, drop = FALSE], con = con)$be ) )
#   suppressWarnings({
#     qa <- nlm(jsreg, ini, y = y, x = x, d = d)
#     qa <- nlm(jsreg, qa$estimate, y = y, x = x, d = d)
#     qa <- nlm(jsreg, qa$estimate, y = y, x = x, d = d)
#   })
#   be <- matrix(qa$estimate, byrow = TRUE, ncol = d)
#   covb <- NULL
#   runtime <- proc.time() - runtime
#   if (B > 1) {
#     betaboot <- matrix( nrow = B, ncol = length(ini) )
#     if (ncores <= 1) {
#       runtime <- proc.time()
#       for (i in 1:B) {
# 	      ida <- rangen::Sample.int(n, n, replace = TRUE)
#         yb <- y[ida, ]
#         xb <- x[ida, ]
#         qa <- nlm(jsreg, be, y = yb, x = xb, d = d)
#         qa <- nlm(jsreg, qa$estimate, y = yb, x = xb, d = d)
#         qa <- nlm(jsreg, qa$estimate, y = yb, x = xb, d = d)
#         betaboot[i, ] <- qa$estimate
#       }  ##  end for (i in 1:B) {
#       covb <- cov(betaboot)
#       runtime <- proc.time() - runtime
#     } else {
#       runtime <- proc.time()
#
#       cl <- parallel::makeCluster(ncores)
#       # Load required packages on workers
#       parallel::clusterEvalQ(cl, library(rangen))
#       # Export only what workers need
#       parallel::clusterExport(cl,
#                              varlist = c("jsreg", "y", "x", "n", "d", "be"),
#                              envir = environment())
#
#       betaboot <- t( parallel::parSapply(cl, 1:B, function(i) {
#         ida <- rangen::Sample.int(n, n, replace = TRUE)
#         yb <- y[ida, ]
#         xb <- x[ida, ]
#         suppressWarnings({
#           qa <- nlm(jsreg, be, y = yb, x = xb, d = d)
#           qa <- nlm(jsreg, qa$estimate, y = yb, x = xb, d = d)
#           qa <- nlm(jsreg, qa$estimate, y = yb, x = xb, d = d)
#           qa$estimate
#         })
#       }))
#
#       parallel::stopCluster(cl)
#       covb <- cov(betaboot)
#       runtime <- proc.time() - runtime
#     }  ##  end (nc <= 1) {
#     nam <- NULL
#     for (i in 1:p)  nam <- c(nam, paste(namy, ":", namx[i], sep = "") )
#     colnames(covb) <- rownames(covb) <- nam
#   }  ##  end if (B > 1) {
#   if ( !is.null(xnew) ) {
#     xnew <- model.matrix(~., data.frame(xnew) )
#     if ( !con )  xnew <- xnew[, -1, drop = FALSE]
#     mu <- cbind( 1, exp(xnew %*% be) )
#     est <- mu / Rfast::rowsums(mu)
#     colnames(est) <- colnames(y)
#   }  else  est <- NULL
#   colnames(be) <- namy
#   rownames(be) <- namx
#   list(runtime = runtime, be = be, covbe = covb, est = est)
# }
