alfa.svm <- function(xnew, y, x, a = seq(-1, 1, by = 0.1), cost = seq(0.2, 2, by = 0.2), gamma = NULL) {

  if ( min(x) == 0 )  a <- a[a > 0]

  if ( is.factor(y) ) {
    task <- "C"
  } else {
    task = "R"
  }

  mod <- list()
  for ( k in 1:length(a) ) {
    z <- Compositional::alfa(x, a[k])$aff
    if ( !is.null(xnew) )  {
      xnew <- as.matrix(xnew)
      if ( dim(xnew)[1] == 1 )  xnew <- t(xnew)
      znew <- alfa(xnew, a[k])$aff
    }
    mod[[ k ]] <- .svmmodel(znew, y, z, task = task, cost = cost, gamma = gamma)
  }

}



.svmmodel <- function(xnew, y, x, task = "R", cost = seq(0.2, 2, by = 0.2), gamma = NULL) {

  if ( is.null(gamma) ) {
    gam <- 1/dim(x)[2]
    gamma <- seq( gam^2, sqrt(gam), length = 10 )
  }
  config <- expand.grid(gamma, cost)
  p <- dim(config)[1]
  mod <- list()

  if ( !is.null(xnew) ) {
    est <- matrix(nrow = dim(xnew)[1], ncol = p)
    xnew <- as.data.frame(xnew)
    colnames(xnew) <- colnames(x)
  }

  if ( task == "R" ) {

    for ( j in 1:p ) {
      mod[[ j ]] <- e1071::svm(y ~., data = as.data.frame(x), type = "eps-regression",
                        gamma = config[j, 1], cost = config[j, 2], scale = FALSE)
      if ( !is.null(xnew) )  est[, j] <- as.numeric( predict(mod[[ j ]], xnew) )
      }

  } else {
    for ( j in 1:p ) {
      mod[[ j ]] <- e1071::svm(y ~., data = as.data.frame(x), type = "C-classification",
                        gamma = config[j, 1], cost = config[j, 2], scale = FALSE)
      if ( !is.null(xnew) )  est[, j] <- as.numeric( predict(mod[[ j ]], xnew) )
    }  ##  end  for ( j in 1:p ) {

  }

  list(mod = mod, config = config, est = est)

}
