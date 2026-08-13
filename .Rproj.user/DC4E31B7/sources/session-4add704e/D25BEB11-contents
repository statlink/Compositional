bic.mixdiri <- function(x, G = 5, tol = 1e-4, graph = FALSE) {
  runtime <- proc.time()
  logn <- log( dim(x)[1] )  ## sample size of the data
  bic <- icl <- 1:G
  mod <- Compositional::diri.nr(x, tol = 1e-6)
  D <- dim(x)[2]
  bic[1] <- icl[1] <-  - 2 * mod$loglik + D * logn  ## BIC assuming one cluster

  for ( vim in 2:G ) {
    a <- Compositional::mix.diri(x, vim, tol = tol)  ## model based clustering for some possible clusters
    d <- dim(a$param)[1]
    nm <- min( table(a$pred) )
    if ( d == vim  &  nm > 10 ) {
      bic[vim] <-  - 2 * a$loglik + ( d - 1 + d * D ) * logn
      icl[vim] <- bic[vim] - 2 * sum( a$probs * log(a$probs), na.rm = TRUE )
    } else  {
      bic[vim] <- bic[vim - 1]
      icl[vim] <- icl[vim - 1]
    }
  }  ## BIC for a range of different clusters

  if ( graph ) {
    ina <- rep(1, G)
    char <- rep(16, G)
    ina[ which.min(icl) ] <- 3  ## chosen number of clusters will appear with red on the plot
    char[ which.min(icl) ] <- 17
    plot(1:G, icl, col = ina, xlab = "Number of components", ylab = "ICL values", cex.lab = 1.3, cex.axis = 1.3)
    abline(v = 1:G, lty = 2, col = "lightgrey")
    abline(h = seq(min(icl, na.rm = FALSE), max(icl, na.rm = FALSE), length = 10), lty = 2, col = "lightgrey" )
    lines(1:G, icl, lwd = 2)
    abline(v = which.min(icl), col = 3, lwd = 2)
    points(1:G, icl, pch = char, col = ina)
  }

  names(bic) <- names(icl) <- paste("g=", 1:G, sep = "")
  runtime <- proc.time() - runtime

  list(bic = bic, icl = icl, runtime = runtime)
}
