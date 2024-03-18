alfalasso.tune <- function(y, x, a = seq(-1, 1, by = 0.1), model = "gaussian", lambda = NULL,
                           type.measure = "mse", nfolds = NULL, folds = NULL) {

  if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds)
  n <- length(y)  ;   f <- numeric(n)
  for ( i in 1:nfolds )   f[ folds[[ i ]] ] <- i

  res <- matrix(nrow = length(a), ncol = 2)
  rownames(a) <- paste("alpha=", a, sep = "")
  colnames(a) <- c("lambda", "performance")

  for ( i in 1:length(a) ) {
    xa <- Compositional::alfa(x, a, h = TRUE)$aff ## apply the alpha-transformation
    mod <- glmnet::cv.glmnet(xa, y, family = model, lambda = lambda, type.measure = type.measure, foldid = f)
    res[i, ] <- c( mod$lambda[ mod$index[1] ], mod$cvm[ mod$index[1] ])
  }
  res
}