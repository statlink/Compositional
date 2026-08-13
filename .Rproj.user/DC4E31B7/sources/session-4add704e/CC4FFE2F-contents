mix.diri <- function(x, g = 2, tol = 1e-4) {

  fun2 <- function(wlika, rswlika, param, x, zx, g, lika) {
    wij <- wlika / rswlika  ## weights
    pj <- Rfast::colmeans(wij)   # PANOS
    for (j in 1:g) {
      param[j, ] <- .wdiri.nr(x, w = wij[, j], zx = zx, tol = 1e-6)
      lika[, j] <- Compositional::ddiri(x, param[j, ]) + log( pj[j] )
    }
    wlika <-  exp(lika) 		#PANOS
    rswlika <- Rfast::rowsums(wlika) #PANOS
    lik <- sum( log( rswlika ) ) 	#PANOS
    wij <- wlika/rswlika
    list(wij = wij, param = param, wlika = wlika, rswlika = rswlika, lika = lika, lik = lik)
  }

  n <- dim(x)[1]  ;  D <- dim(x)[2]
  zx <- t( Rfast::Log(x) )

  lik <- NULL
  lika <- matrix(nrow = n, ncol = g)
  param <- matrix(nrow = g, ncol = D)

  runtime <- proc.time()
  ## Step 1

  if ( g > 1 ) {
    cl <- kmeans(x, g, nstart = 20)$cl
  } else  cl <- rep(1, n)
  wij <- tabulate(cl)

  while ( min(wij) <= 10 ) {
    g <- g - 1
    lika <- matrix(nrow = n, ncol = g)
    param <- matrix(nrow = g, ncol = 5)
    cl <- kmeans(x, g, nstart = 20)$cl
    wij <- tabulate(cl)
  }

  for ( j in 1:g ) {
    mod <- Compositional::diri.nr(x[cl == j, ], tol = 1e-6)
    param[j, ] <-  mod$param
    lika[, j] <- Compositional::ddiri(x, param[j, ])
  }

  wlika <- exp(lika)
  rswlika <- Rfast::rowsums(wlika)

  ep <- fun2(wlika, rswlika, param, x, zx, g, lika)
  lik[1] <- ep$lik
  ep2 <- fun2(ep$wlika, ep$rswlika, ep$param, x, zx, g, ep$lika)
  lik[2] <- ep2$lik

  i <- 2
  while ( lik[i] - lik[i - 1] > tol ) {
    i <- i + 1
    ep <- ep2
    ep2 <- fun2(ep$wlika, ep$rswlika, ep$param, x, zx, g, ep$lika)
    lik[i] <- ep2$lik
  }
  res <- ep2
  if ( ep$lik > ep2$lik )  res <- ep

  pj <- Rfast::colmeans(res$wij)
  loglik <- res$lik   #sum( log( Rfast::colsums( pj * t( exp( res$lika ) ) ) ) )
  ta <- Rfast::rowMaxs(res$wij)  ## estimated cluster of each observation
  param <- cbind(pj, res$param)
  runtime <- proc.time() - runtime
  colnames(param) <- c( "probs", paste("alpha", 1:D, sep = "") )
  rownames(param) <- paste("Cluster", 1:g, sep = " ")

  list( param = param, loglik = loglik, probs = res$wij, pred = ta, iters = i, runtime = runtime )
}



.wdiri.nr <- function (x, w, zx, tol = 1e-06) {
  dm <- dim(x)
  n <- dm[1]
  p <- dm[2]
  W <- sum(w)
  m <- Rfast::eachcol.apply(x, w) / W
  gm <- as.vector(zx %*% w)
  down <-  - sum( m * ( gm / W - log(m) ) )
  sa <- 0.5 * (p - 1) / down
  a1 <- sa * m
  z <- W * digamma(sa)
  g <- z - W * digamma(a1) + gm
  qk <-  - W * trigamma(a1)
  b <- sum(g / qk) / ( 1 / z - sum(1 / qk) )
  a2 <- a1 - (g - b) / qk
  while ( sum( abs(a2 - a1 ) ) > tol ) {
    a1 <- a2
    z <- W * digamma( sum(a1) )
    g <- z - W * digamma(a1) + gm
    qk <-  - W * trigamma(a1)
    b <- sum(g / qk) / ( 1 / z - sum(1 / qk) )
    a2 <- a1 - (g - b) / qk
  }
  a2
}
