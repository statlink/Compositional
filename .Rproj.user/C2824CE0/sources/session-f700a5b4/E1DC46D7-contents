ascrq <- function(y, x, xnew = NULL) {

  py <- dim(y)[2]   ;    px <- dim(x)[2]
  pyx <- py * px    ;    n <- dim(y)[1]
  n <- dim(y)[1]    ;    npy <- n * py

  X <- matrix(0, npy, pyx)
  indr <- matrix( 1:npy, ncol = py )
  indc <- matrix( 1:pyx, ncol = py )
  for ( i in 1:py )  X[ indr[, i], indc[, i] ] <- x
  Y <- as.vector(y)

  R <- NULL
  for (i in 1:py)  R <- cbind(R, diag(px))
  R <- rbind( R, -R, diag(pyx), -diag(pyx) )
  r <- c( rep(1, px), rep(-1, px), rep(0, pyx), rep(-1, pyx) )

  res <- list()
  for ( j in 1:length(a) ) {
    ya <- y^a[j]
    ya <- ya / Rfast::rowsums(ya)
    mod <- quantreg::rq(Y ~ X - 1, data = data.frame(Y = Y, X = X), method = "fnc", R = R, r = r)
    be <- matrix( coef(a), ncol = py)
    est <- ( xnew %*% be )^( 1/a[j] )
    est <- est / Rfast::rowsums(est)
    res[[ j ]] <- est
  }

  res
}
