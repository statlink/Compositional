tflr.indeptest <- function(y, x, tol = 1e-6, R = 999, ncores = 1) {
  kl <- Compositional::tflr.irls(y, x)$kl
  n <- dim(y)[1]
  pkl <- numeric(R)
  
  if ( ncores <= 1 ) {
    for ( i in 1:R ) {
      id <- Rfast2::Sample.int(n, n)
      pkl[i] <- Compositional::tflr.irls(y, x[id, ], tol = tol)$kl
    }
  } else {
    cl <- parallel::makeCluster(ncores)
    # Load required packages on all workers
    parallel::clusterEvalQ(cl, {
      library(Rfast2)
      library(Compositional)
    })
    # Export only what workers need
    parallel::clusterExport(cl, varlist = c("y", "x", "n", "tol"), 
                           envir = environment())
    
    pkl <- parallel::parSapply(cl, 1:R, function(i) {
      id <- Rfast2::Sample.int(n, n)
      Compositional::tflr.irls(y, x[id, ], tol = tol)$kl
    })
    
    parallel::stopCluster(cl)
  }
  
  ( sum(pkl < kl) + 1 ) / (R + 1)
}







