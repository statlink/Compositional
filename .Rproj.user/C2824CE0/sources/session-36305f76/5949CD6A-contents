eqdist.etest <- function(x, y, R = 999) {
  nx <- dim(x)[1]  ;  ny <- dim(y)[1]
  n <- nx + ny
  stat <- Rfast::Dist(x, result = "sum") / nx + Rfast::Dist(y, result = "sum") / ny
  z <- rbind(x, y) 
  pstat <- numeric(R)
  for ( i in 1:R ) {
    id <- Rfast2::Sample.int(n, nx)
    pstat[i] <- Rfast::Dist(z[id, ], result = "sum") / nx + 
             Rfast::Dist(z[-id, ], result = "sum") / ny
  }
  ( sum( pstat >= stat ) + 1 ) / (R + 1)
}     