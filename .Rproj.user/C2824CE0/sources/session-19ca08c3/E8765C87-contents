alfasvm.tune(y, x, a, type, cost = seq(0.2, 2, by = 0.2), gamma = NULL, ncores = 1,
             folds = NULL, nfolds = 10, stratified = TRUE, seed = NULL)




.svm.tune <- function(y, x, task = "R", cost = seq(0.2, 2, by = 0.2), gamma = NULL,
                      folds = NULL, nfolds = 10, stratified = TRUE, seed = NULL) {

  if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds, statified = stratified, seed = seed )
  nfolds <- length(folds)
  if ( is.null(gamma) ) {
    gam <- 1/dim(x)[2]
    gamma <- seq( gam^2, sqrt(gam), length = 10 )
  }
  config <- expand.grid(gamma, cost)
  p <- dim(config)[1]
  per <- matrix(nrow = nfolds, ncol = p)
  if ( is.factor(y) ) {
    type <- "C"
  } else type = "R"

  if ( type == "R" ) {

    runtime <- proc.time()
    for ( k in 1:10 ) {
      ytrain <- y[ -folds[[ k ]] ]
      ytest <- y[ folds[[ k ]] ]
      xtrain <- x[-folds[[ k ]], ]
      xtest <- as.data.frame( x[folds[[ k ]], ] )
      colnames(xnew) <- colnames(xtrain)
      st <- matrix(nrow = length(ytest), ncol = p)

      for ( j in 1:p ) {
        mod <- e1071::svm(ytrain ~., data = as.data.frame(xtrain), type = "eps-regression",
                         gamma = config[j, 1], cost = config[j, 2], scale = FALSE)
        st[, j] <- as.numeric( predict(mod, xnew) )
      }  ##  end  for ( j in 1:p ) {
      per[k, ] <- Rfast2::colmses(ytest, st)

    }  ##  end  for (k in 1:nfolds) {

    runtime <- proc.time() - runtime
    per <- cbind(config, Rfast::colmeans(per) )
    colnames(per) <- c("gamma", "cost", "mse")
    ind <- which.min(per[, 3])

  } else {

    for ( k in 1:10 ) {
      ytrain <- y[ -folds[[ k ]] ]
      ytest <- y[ folds[[ k ]] ]
      xtrain <- x[-folds[[ k ]], ]
      xtest <- as.data.frame( x[folds[[ k ]], ] )
      colnames(xnew) <- colnames(xtrain)
      st <- matrix(nrow = length(ytest), ncol = p)

      for ( j in 1:p ) {
        mod <- e1071::svm(ytrain ~., data = as.data.frame(xtrain), type = "C-classification",
                          gamma = config[j, 1], cost = config[j, 2], scale = FALSE)
        st[, j] <- as.numeric( predict(mod, xnew) ) - 1
      }  ##  end  for ( j in 1:p ) {
      per[k, ] <- Rfast::colaccs(ytest, st)

    }  ##  end  for (k in 1:nfolds) {

    per <- cbind(config, Rfast::colmeans(per) )
    colnames(per) <- c("gamma", "cost", "mse")
    ind <- which.max(per[, 3])

  }

  list(per = per, perf = per[ind, ])
}
