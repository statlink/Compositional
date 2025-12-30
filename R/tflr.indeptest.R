tflr.indeptest <- function(y, x, tol = 1e-6, R = 999, ncores = 1) {

  kl <- Compositional::tflr(y, x)$kl
  n <- dim(y)[1]
  pkl <- numeric(R)

  if ( ncores <= 1 ) {
    for ( i in 1:R ) {
      id <- Rfast2::Sample.int(n, n)
      pkl[i] <- Compositional::tflr(y, x[id, ], tol = tol)$kl
    }

  } else {
    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport( cl, c("y", "x", "n"), envir = environment() )
    pkl <- parallel::parSapply(cl, 1:R, function(i) {
      id <- Rfast2::Sample.int(n, n)
      Compositional::tflr(y, x[id, ], tol = tol)$kl
    })
    parallel::stopCluster(cl)    
  }

 ( sum(pkl < kl) + 1 ) / (R + 1)
}










