eqdist.etest <- function(x, y, R = 999) {
  nx <- dim(x)[1]  ;  ny <- dim(y)[1]
  n <- nx + ny
  z <- rbind(x, y)
  
  tot <- Rfast::Dist(z, result = "sum")
  mij <- Rfast::dista(x, y, result = "sum")
  if ( nx > ny ) { 
    dy <- Rfast::Dist(y, result = "sum") 
    dx <- tot - mij - dy 
  } else {
    dx <- Rfast::Dist(x, result = "sum") 
    dy <- tot - mij - dx 
  }
  stat <- mij - ny * dx/nx -  nx * dy/ny

  pstat <- numeric(R)
  for ( i in 1:R ) {
    id <- Rfast2::Sample.int(n, nx)
    zx <- z[id, ]  ;  zy <- z[-id, ]
    mij <- Rfast::dista(zx, zy, result = "sum")
    if ( nx > ny ) { 
      dy <- Rfast::Dist(zy, result = "sum") 
      dx <- tot - mij - dy 
    } else {
      dx <- Rfast::Dist(zx, result = "sum") 
      dy <- tot - mij - dx 
    }
    pstat[i] <- mij - ny * dx/nx -  nx * dy/ny 
  }
  ( sum( pstat >= stat ) + 1 ) / (R + 1)
}




# eqdist.etest <- function(x, y, R = 999) {
  # nx <- dim(x)[1]  ;  ny <- dim(y)[1]
  # n <- nx + ny
  # stat <- Rfast::edist(x, y)
  # z <- rbind(x, y)
  # pstat <- numeric(R)
  # for ( i in 1:R ) {
    # id <- Rfast2::Sample.int(n, nx)
    # pstat[i] <- Rfast::edist(z[id, ], z[-id, ])
  # }
  # ( sum( pstat >= stat ) + 1 ) / (R + 1)
# }
