beta.reg.irls <- function(y, x, xnew = NULL, tol = 1e-6, maxit = 100) {
  x <- model.matrix(y~., data = as.data.frame(x) ) 
  n <- dim(x)[1]   ;   p <- dim(x)[2]

  ystar  <- qlogis(y)
  ilogit <- function(eta) 1 / ( 1 + exp(-eta) )
  ly <- log(y)   ;   ly1 <- log(1 - y)
  beta <- solve( crossprod(x), crossprod(x, y) )
  mu <- drop( ilogit(x %*% beta) )          # <- as.vector here
  phi <- max(  mean(mu * (1 - mu) ) / var(y) - 1, 2 )
  lphi <- log(phi)

  loglik.old <-  -Inf

  for ( it in 1:maxit ) {
    mu  <- drop( ilogit(x %*% beta) )          # <- and here
    phi <- exp(lphi)
    mustar <- digamma(mu * phi) - digamma( (1 - mu) * phi )
    rphi_i <- mu * (ystar - mustar) + ly1 - digamma( (1 - mu) * phi ) + digamma(phi)
    dmu.deta <- mu * (1 - mu)
    Ubeta <- drop( crossprod(x, dmu.deta * phi * (ystar - mustar)) )
    Ulphi <- phi * sum(rphi_i)
    psi1_mp  <- trigamma(mu * phi)
    psi1_1mp <- trigamma( (1 - mu) * phi )
    psi1_phi <- trigamma(phi)
    Wbb_i <- phi^2 * (psi1_mp + psi1_1mp) * dmu.deta^2
    Wbg_i <- phi * ( psi1_mp * mu - psi1_1mp * (1 - mu) ) * dmu.deta * phi
    Wgg <- phi^2 * sum( psi1_mp * mu^2 + psi1_1mp * (1 - mu)^2 - psi1_phi )

    Ibb <- crossprod(x, x * Wbb_i)     # now fine: Wbb_i is a plain vector
    Ibg <- drop( crossprod(x, Wbg_i) )
    I <- rbind( cbind(Ibb, Ibg), c(Ibg, Wgg) )
    U <- c(Ubeta, Ulphi)
    step <- solve(I, U)
    beta <- beta + step[1:p]
    lphi <- lphi + step[p+1]

    mu <- drop( ilogit(x %*% beta) )
    phi <- exp(lphi)
    loglik <- sum( lgamma(phi) - lgamma(mu * phi) - lgamma( (1-mu) * phi ) +
                   (mu * phi - 1) * ly + ( (1 - mu) * phi - 1 ) * ly1 )

    if ( abs(loglik - loglik.old) < tol )  break
    loglik.old <- loglik
  }
   
  s <- solve(I)
  sephi <- a$phi * sqrt( s[p + 1, p + 1] )

  be <- beta
  se <- sqrt( diag(s)[1:p] ) 
  stat <- (be/se)^2
  pvalue <- pchisq(stat, 1, lower.tail = FALSE)
  info <- cbind(be, se, stat, pvalue)
  rownames(info) <- colnames(x)

  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )
    est <- exp( - as.vector( xnew %*% be[, 1] ) )
    est <- 1 / (1 + est)
  } else  est <- NULL


  list( iters = i, loglik = loglik, info = info, phi = phi, sephi = sephi, est = est )
}