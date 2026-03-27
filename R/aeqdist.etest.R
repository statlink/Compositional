aeqdist.etest <- function(x, sizes, a = 1, R = 999) {

  x1 <- x[1:sizes[1], ]
  x2 <- x[-c(1:sizes[1]), ]

  if ( length(a) == 1 ) {
    y1 <- Rfast::standardise( Compositional::alfa(x1, a)$aff )
    y2 <- Rfast::standardise( Compositional::alfa(x2, a)$aff )
    res <- Compositional::eqdist.etest(y1, y2, R = R)
  } else {
    if ( min(x) == 0 )  a <- a[ a > 0 ]
    len <- length(a)
    res <- numeric(len)
    for ( i in 1:len )  {
      y1 <- Rfast::standardise( Compositional::alfa(x1, a[i])$aff )
      y2 <- Rfast::standardise( Compositional::alfa(x2, a[i])$aff )
      res[i] <- Compositional::eqdist.etest(y1, y2, R = R)
    }
    names(res) <- a
  }

  res
}
