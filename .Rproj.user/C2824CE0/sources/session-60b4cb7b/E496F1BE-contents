alfasvm.tune <- function(y, x, a = seq(-1, 1, by = 0.1), cost = seq(0.2, 2, by = 0.2), gamma = NULL, ncores = 1,
                         folds = NULL, nfolds = 10, stratified = TRUE, seed = NULL, graph = FALSE) {

  if ( min(x) == 0 )  a <- a[a > 0]

  if (ncores > 1) {
    runtime <- proc.time()
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    if ( is.factor(y) ) {
      task <- "C"
      if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds, stratified = stratified, seed = seed )
      nfolds <- length(folds)
    } else {
      task <- "R"
      if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds, stratified = FALSE, seed = seed )
      nfolds <- length(folds)
    }
    ww <- foreach(k = 1:length(a), .combine = cbind, .export = c(".svmtune", "svm", "colaccs",
                  "colmses", "colmeans"), .packages = c("e1071", "Rfast", "Rfast2") ) %dopar% {
     z <- Compositional::alfa(x, a[k])$aff
     per <- as.numeric( .svmtune(y, z, task = task, cost = cost, gamma = gamma, folds = folds)$perf )
     return(per)
   }
   parallel::stopCluster(cl)
   per <- ww
   runtime <- proc.time() - runtime

  } else {
    if ( is.factor(y) ) {
      task <- "C"
      if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds, stratified = stratified, seed = seed)
      nfolds <- length(folds)
    } else {
      task <- "R"
      if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds, stratified = FALSE, seed = seed)
      nfolds <- length(folds)
    }
    per <- matrix(nrow = length(a), ncol = 3)
    runtime <- proc.time()

    for ( k in 1:length(a) ) {
      z <- Compositional::alfa(x, a[k])$aff
      per[k, ] <- as.numeric( .svmtune(y, z, task = task, cost = cost, gamma = gamma, folds = folds)$perf )
    }
    runtime <- proc.time() - runtime
  }

  if (graph) {
    plot(a, per[, 3], type = "b", ylim = c( min(per[, 3]), max(per[, 3]) ), ylab = "Estimated performance",
         xlab = expression( paste(alpha, " values") ), cex.lab = 1.2, cex.axis = 1.2, pch = 16, col = "green")
    abline(v = a, col = "lightgrey", lty = 2)
    abline(h = seq(min(per[ ,3]), max(per[, 3]), length = 10), col = "lightgrey", lty = 2)
  }

  rownames(per) <- paste("alpha=", a, sep = "")
  colnames(per) <- c("gamma", "cost", "performance")

  if (task == "C") {
    ind <- which.max(per[, 3])
  } else ind <- which.min(per[, 3])

  list(per = per, performance = per[ind, 3], best_a = a[ind], runtime = runtime)
}





.svmtune <- function(y, x, task = "R", cost = seq(0.2, 2, by = 0.2), gamma = NULL,
                      folds = NULL) {

  nfolds <- length(folds)
  if ( is.null(gamma) ) {
    gam <- 1/dim(x)[2]
    gamma <- seq( gam^2, sqrt(gam), length = 10 )
  }
  config <- expand.grid(gamma, cost)
  p <- dim(config)[1]
  per <- matrix(nrow = nfolds, ncol = p)

  if ( task == "R" ) {

    runtime <- proc.time()
    for ( k in 1:nfolds ) {
      ytrain <- y[ -folds[[ k ]] ]
      ytest <- y[ folds[[ k ]] ]
      xtrain <- x[-folds[[ k ]], ]
      xtest <- as.data.frame( x[folds[[ k ]], ] )
      colnames(xtest) <- colnames(xtrain)
      st <- matrix(nrow = length(ytest), ncol = p)

      for ( j in 1:p ) {
        mod <- e1071::svm(ytrain ~., data = as.data.frame(xtrain), type = "eps-regression",
                         gamma = config[j, 1], cost = config[j, 2], scale = FALSE)
        st[, j] <- as.numeric( predict(mod, xtest) )
      }  ##  end  for ( j in 1:p ) {
      per[k, ] <- Rfast2::colmses(ytest, st)

    }  ##  end  for (k in 1:nfolds) {

    runtime <- proc.time() - runtime
    per <- cbind(config, Rfast::colmeans(per) )
    colnames(per) <- c("gamma", "cost", "mse")
    ind <- which.min(per[, 3])

  } else {

    for ( k in 1:nfolds ) {
      ytrain <- y[ -folds[[ k ]] ]
      ytest <- y[ folds[[ k ]] ]
      xtrain <- x[-folds[[ k ]], ]
      xtest <- as.data.frame( x[folds[[ k ]], ] )
      colnames(xtest) <- colnames(xtrain)
      st <- matrix(nrow = length(ytest), ncol = p)

      for ( j in 1:p ) {
        mod <- e1071::svm(ytrain ~., data = as.data.frame(xtrain), type = "C-classification",
                          gamma = config[j, 1], cost = config[j, 2], scale = FALSE)
        st[, j] <- as.numeric( predict(mod, xtest) )
      }  ##  end  for ( j in 1:p ) {
      per[k, ] <- Rfast2::colaccs( as.numeric(ytest), st)

    }  ##  end  for (k in 1:nfolds) {

    per <- cbind(config, Rfast::colmeans(per) )
    colnames(per) <- c("gamma", "cost", "acc")
    ind <- which.max(per[, 3])

  }

  list(per = per, perf = per[ind, ])
}
