alfa.tune <- function(x, B = 1, ncores = 1) {
  n <- dim(x)[1]  ## sample size
  fc <- (n - 1) / n
  D <- dim(x)[2]   ;   d <- D - 1
  con <-  - 0.5 * n * d * log(2 * pi) - 0.5 * n * d

  pa <- function(a, x, n, fc) {
    trans <- Compositional::alfa(x, a)
    z <- trans$aff  ## the alpha-transformation
    #- 0.5 * n * log( abs( det( fc * cov(z) ) ) ) + trans$sa
    - 0.5 * n * as.numeric( determinant( fc * cov(z), logarithm = TRUE)$modulus ) + trans$sa
  }

  if ( B == 1 ) {
    suppressWarnings({
      ell <- optimize(pa, c(-1, 1), x = x, n = n, fc = fc, maximum = TRUE )
    })
    aff0 <- Compositional::alfa(x, 0)
    z0 <- aff0$aff
    lik0 <-  - 0.5 * n * log( abs( det( fc * Rfast::cova(z0) ) ) ) + aff0$sa
    result <- c(ell$maximum, ell$objective + con, lik0 + con)
    names(result) <- c("best alpha", "max log-lik", "log-lik at 0")

  } else {  ## bootstrap confidence intervals
    suppressWarnings({
      ell <- optimize(pa, c(-1, 1), x = x, n = n, fc = fc, maximum = TRUE )
    })
    ab <- numeric(B)

    if ( ncores <= 1 ) {
      runtime <- proc.time()
      for (i in 1:B) {
        ind <- rangen::Sample.int(n, n, replace = TRUE)
        suppressWarnings({
          ab[i] <- optimize(pa, c(-1, 1), x = x[ind, ], n = n, fc = fc, maximum = TRUE )$maximum
        })
      }
      runtime <- proc.time() - runtime

    } else {
      runtime <- proc.time()
      cl <- parallel::makeCluster(ncores)
      # Load required packages on workers
      parallel::clusterEvalQ(cl, {
        library(rangen)
        library(Compositional)
      })
      # Export only what workers need
      parallel::clusterExport(cl,
                             varlist = c("pa", "x", "n", "D", "fc"),
                             envir = environment())

      ab <- parallel::parSapply(cl, 1:B, function(i) {
        ind <- rangen::Sample.int(n, n, replace = TRUE)
        suppressWarnings({
          optimize(pa, c(-1, 1), x = x[ind, ], n = n, fc = fc, maximum = TRUE )$maximum
        })
      })

      parallel::stopCluster(cl)
      runtime <- proc.time() - runtime
    }

    param <- c(ell$maximum, ell$objective + con, quantile( ab, c(0.025, 0.975) ) )
    names(param)[1:2] <- c("best alpha", "max log-lik")
    hist( ab, main = "Bootstrapped alpha values", xlab = expression( paste(alpha, " values", sep = "") ),
         cex.lab = 1.2, cex.axis = 1.2 )
    abline(v = ell$maximum, col = 3)
    abline(v = mean(ab), lty = 2, col = 4)
    message <- paste("The green is the best alpha value. The blue line is the bootstrap mean value of alpha.")
    result <- list(param = param, message = message, runtime = runtime )
  }  ##  end if ( B == 1 ) {

  result
}
