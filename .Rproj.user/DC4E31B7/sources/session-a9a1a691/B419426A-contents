alfa.profile <- function(x, a = seq(-1, 1, by = 0.01) ) {
  D <- dim(x)[2]   ;   d <- D - 1
  n <- dim(x)[1]
  f <- (n - 1) / n
  qa <- numeric( length(a) )
  con <-  - 0.5 * n * d * log(2 * pi) - 0.5 * n * d

  for ( i in 1:length(a) ) {
    trans <- Compositional::alfa(x, a[i])
    aff <- trans$aff
    qa[i] <-  - 0.5 * n * as.numeric( determinant( f * cov(aff), logarithm = TRUE)$modulus ) + trans$sa
  }
  qa <- qa + con

  ## the green lines show a 95% CI for the true value of
  ## alpha using a chi-square distribution
  b <- max(qa) - qchisq(0.95, 1)/2
  plot(a, qa, type = "l", xlab = expression( paste(alpha, " values", sep = "") ),
  ylab = "Profile log-likelihood", cex.axis = 1.2, cex.lab = 1.2)
  abline(h = b, col = 2)
  ci <- c( min(a[qa >= b]), max(a[qa >= b]) )
  names(ci) <- paste(c("2.5", "97.5"), "%", sep = "")
  abline(v = ci[1], col = 3, lty = 2)
  abline(v = ci[2], col = 3, lty = 2)
  res <- c(a[which.max(qa)], max(qa), qa[a == 0])
  names(res) <- c('alfa', 'max.log.lik', 'log.lik0')
  list(result = res, ci = ci)
}
