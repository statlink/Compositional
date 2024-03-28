tflr.indeptest <- function(y, x, B = 999, ncores = 1) {

  kl <- tflr2(y, x)$kl
  n <- dim(y)[1]
  pkl <- numeric(B)

  if (ncores <= 1) {
    for (i in 1:R) {
      id <- Rfast2::Sample.int(n, n)
      pkl[i] <- Compositional::tflr(y, x[id, ])$kl
    }
	
  } else {  
    requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    pkl <- foreach::foreach(i = 1:R, .combine = "c", 
                   .packages = c("Compositional", "Rfast", "Rfast2") %dopar% {
      id <- Rfast2::Sample.int(n, n)
      return( Compositional::tflr(y, x[id, ])$kl )
    }
  }

  ( sum(pkl < kl) + 1 ) / (B + 1)
}

