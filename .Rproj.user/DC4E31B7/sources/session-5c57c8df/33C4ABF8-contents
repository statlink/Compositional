dzad <- function(y, phi, mu, logged = TRUE) {

  dm <- dim(y)
  n <- dm[1]  ;  D <- dm[2]  ;  d <- D - 1
  f <- numeric(n)
  a1 <- which( Rfast::rowsums( y > 0 ) == D )
  a2 <- which( Rfast::rowsums( y > 0 ) != D )
  n1 <- length(a1)
  n2 <- n - n1

  za <- y[a2, , drop = FALSE]
  za[za == 0] <- 1
  za[ za < 1 ] <- 0
  theta <- apply(za, 1, paste, collapse = ",")
  freq <- table(theta)
  p2 <- as.numeric(freq[theta]) / n
  p1 <- n1 / n

  y1 <- y[a1, , drop = FALSE]
  ly1 <- log(y1)
  ly2 <- log( y[a2, , drop = FALSE] )

  mu <- matrix(mu, nrow = n1, ncol = D, byrow = TRUE)
  mu2 <- mu[1:n2, ]
  ly3 <- ly2
  ind <- which(is.infinite(ly2))
  ly3[ind] <- 0
  mu2[ind] <- 0

  w <- lgamma( phi * mu2 )
  w[is.infinite(w)] <- 0
  f[a2] <- lgamma(phi) - Rfast::rowsums(w) + Rfast::rowsums( (mu2 * phi - 1) * ly3 ) + log(p2)
  ba <- phi * mu
  f[a1] <- lgamma(phi) - sum( lgamma( ba[1, ] ) ) + Rfast::rowsums( (ba - 1) * ly1 ) + log(p1)

  if ( !logged )
  f <- exp(f)
  f

}


