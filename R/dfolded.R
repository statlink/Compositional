dfolded <- function(x, a, p, mu, su, logged = TRUE) {

  if ( is.null(dim(x)[1]) )  x <- matrix(x, nrow = 1)
  d <- dim(x)[2] - 1
  h <- t( helm(d + 1) )
  down <- 1/sqrt( det( 2 * pi * su) )

  z1 <- Compositional::alef(x, a)$aff
  y1 <- z1 %*% h
  if ( a > 0 ) {
    lam <- 1 / ( a^2 * Rfast::rowMins(y, value = TRUE)^2 )    ##  1/apply(a * y, 1, min)^2
  } else {
    lam <- 1 / ( a^2 * Rfast::rowMaxs(y, value = TRUE)^2 )    ##  1/apply(a * y, 1, max)^2
  }

  y2 <- lam * y1
  f <- p * down * exp( -0.5 * Rfast::mahala(y1, mu, su) ) + (1 - p) * down * lam^d * exp( -0.5 * Rfast::mahala(y2, mu, su) )
  if (logged) f <- log(f)
  f
}
