kernreg.tune <- function(y, x, h = seq(0.1, 1, length = 10), type = "gauss",
                         nfolds = 10, folds = NULL, seed = NULL, graph = FALSE, ncores = 1) {

  if ( is.matrix(y) )  {
    n <- dim(y)[1]
    D <- dim(y)[2]
  } else {
    n <- length(y)
    D <- 1
  }

  runtime <- proc.time()
  if ( !is.matrix(x) )  x <- as.matrix(x)
  ina <- 1:n
  if ( is.null(folds) )  folds <- Compositional::makefolds(ina, nfolds = nfolds,
                                  stratified = FALSE, seed = seed)
  nfolds <- length(folds)
  msp <- matrix( nrow = length(folds), ncol = length(h) )

  if ( ncores <= 1 ) {
    for ( i in 1:nfolds ) {
      nu <- folds[[ i ]]
      if ( D > 1 )  {
        ytrain <- y[-nu, , drop = FALSE]
        ytest <- y[nu, , drop = FALSE]
      } else {
        ytest <- y[nu]
        ytrain <- y[-nu]
      }
      xtest <- x[nu, , drop = FALSE]
      xtrain <- x[-nu, , drop = FALSE]
      est <- Compositional::kern.reg(xtest, ytrain, xtrain, h, type = type)
      if ( D > 1 ) {
        for ( j in 1:length(h) )  msp[i, j] <- mean( est[[ j ]] - ytest)^2
      } else  msp[i, ] <- Rfast::colmeans( (ytest - est)^2 )
    }

  } else {
    pe <- numeric( length(h) )
    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport( cl, varlist = ls(), envir = environment() )
    msp <- t( parallel::parSapply(cl, 1:nfolds, function(i) {
      nu <- folds[[ i ]]
      if ( D > 1 )  {
        ytrain <- y[-nu, , drop = FALSE]
        ytest <- y[nu, , drop = FALSE]
      } else {
        ytest <- y[nu]
        ytrain <- y[-nu]
      }
      xtest <- x[nu, , drop = FALSE]
      xtrain <- x[-nu, , drop = FALSE]
      est <- Compositional::kern.reg(xtest, ytrain, xtrain, h, type = type)
      if ( D > 1 ) {
        for ( j in 1:length(h) )  pe[j] <- mean( est[[ j ]] - ytest)^2
      } else  pe <- Rfast::colmeans( (ytest - est)^2 )
      pe
    }))
    parallel::stopCluster(cl)
  }

  runtime <- proc.time() - runtime

  mspe <- Rfast::colmeans(msp)
  names(mspe) <- paste("h=", h, sep = "")
  if (graph) {
    plot(h, mspe, xlab = "Bandwidth parameter (h)", ylab = "MSPE", cex.lab = 1.2, cex.axis = 1.2)
    abline(v = h, col = "lightgrey", lty = 2)
    abline(h = seq(min(mspe), max(mspe), length = 10), col = "lightgrey", lty = 2)
    points(h, mspe, pch = 19)
    lines(h, mspe, lwd = 2)
  }

  list(mspe = mspe, h = h[ which.min(mspe) ], performance = min(mspe), runtime = runtime)
}