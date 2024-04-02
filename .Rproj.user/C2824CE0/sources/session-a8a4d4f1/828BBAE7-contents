tflr.indeptest <- function(y, x, B = 999) {

  kl <- tflr2(y, x)$kl
  n <- dim(y)[1]
  pkl <- numeric(B)
  for (i in 1:R) {
    id <- Rfast2::Sample.int(n, n)
    pkl[i] <- Compositional::tflr(y, x[id, ])$kl
  }

  ( sum(pkl < kl) + 1 ) / (B + 1)
}

