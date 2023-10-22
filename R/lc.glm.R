lc <- function(y, x, z = NULL, xnew = NULL, znew = NULL) {

  p <- dim(x)[2]

  poisdev <- function(be) {
    est <- x %*% be
     -2 * sum( y * est, na.rm = TRUE )
  }

    beini <- Compositional::lc.reg(y, x)$be
    x <- cbind( 1, log(x) )

    con <- function(be){
      f <- sum(be[-1])
      list(ceq = f, c = NULL)
    }

    dev <- poisdev

    runtime <- proc.time()
    f1 <- NlcOptim::solnl( X = beini, dev, con, lb = rep(-200, p + 1), ub = rep(200, p + 1) )
    f2 <- NlcOptim::solnl( f1$par, dev, con, lb = rep(-200, p + 1), ub = rep(200, p + 1) )
    while ( f1$fn - f2$fn > 1e-04 ) {
      f1 <- f2
      f1 <- NlcOptim::solnl( f2$par, dev, con, lb = rep(-200, p + 1), ub = rep(200, p + 1) )
      f2 <- NlcOptim::solnl( f1$par, dev, con, lb = rep(-200, p + 1), ub = rep(200, p + 1) )
    }
    runtime <- proc.time() - runtime

    devi <- f2$fn + 2 * sum( y * log(y), na.rm = TRUE ) 
    be <- as.vector(f2$par)
    if ( is.null( colnames(x) ) ) {
      names(be) <- c( "constant", paste("X", 1:p, sep = "") )
    } else  names(be) <- colnames(x)

    est <- NULL
    if ( !is.null(xnew) ) {
      est <- cbind(1, log(xnew) ) %*% be
      est <- exp(est)
    }

    res <- list( runtime = runtime, devi = devi, be = be, est = est )

}
