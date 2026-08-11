zadr.irls <- function(y, x, xnew = NULL, tol = 1e-6, maxit = 100) {
  ## y is the compositional data (can contain structural zeros)
  ## x is the independent variable(s)
  ##
  ## For a row with no zeros the usual D-dimensional Dirichlet log-density
  ## is used. For a row with zeros in categories Z (survivors P = {1..D}\Z),
  ## zadr()'s .mixreg renormalises the softmax mean over the SURVIVORS only
  ## (mu2 <- exp(x2 %*% be); mu2[zero-cols] <- 0; mu2 <- mu2/rowSums(mu2))
  ## and fits a genuine |P|-dimensional Dirichlet on the survivors with the
  ## same total concentration phi. This function reproduces that exactly,
  ## via Fisher scoring instead of optim.

  n <- dim(y)[1]   ;  D <- dim(y)[2]  ;  p <- D - 1
  x <- model.matrix( y~., data = as.data.frame(x) )
  K <- dim(x)[2]

  ## -- multinomial "zero-pattern" constant (independent of beta/phi,
  ##    added at the end so loglik matches zadr()'s scale) --
  pos <- y > 0
  a1 <- which( Rfast::rowsums( pos ) == D )
  a2 <- which( Rfast::rowsums( pos ) != D )
  n1 <- length(a1)   ;   n2 <- n - n1
  if ( n2 > 0 ) {
    za <- y[a2, , drop = FALSE]
    za[za == 0] <- 1  ;  za[za < 1] <- 0
    theta <- as.vector( table( apply(za, 1, paste, collapse = ",") ) )
    const <- n1 * log(n1 / n) + sum( theta * log(theta / n) )
  } else  const <- n1 * log(n1 / n)

  ## survivor index set per row
  Pidx <- lapply( 1:n, function(i) which( pos[i, ] ) )
  logy <- log(y)

  beta <- matrix(0, nrow = K, ncol = p)   # columns = components 2..D
  phi <- 1.0
  loglik_old <-  -Inf

  txi <- list()
  for ( i in 1:n ) txi[[ i ]] <- tcrossprod(x[i, ])

  for ( iter in 1:maxit ) {
    eta <- x %*% beta                     # n x p
    exp_eta <- exp(eta)
    mu <- cbind(1, exp_eta) / (1 + Rfast::rowsums(exp_eta))   # n x D

    psi1_phi <- trigamma(phi)
    loglik <- 0
    S_vec <- numeric(K * p)
    I_mat <- matrix(0, K * p, K * p)
    S_phi <- 0
    H_phi <- 0

    for ( i in 1:n ) {
      Pi <- Pidx[[ i ]]
      x_i <- x[i, ]
      mu_i <- mu[i, ]
      Mi <- sum( mu_i[Pi] )
      mu2 <- mu_i[Pi] / Mi            # renormalised over survivors
      alpha2 <- phi * mu2
      y_log <- logy[i, Pi]

      loglik <- loglik + lgamma(phi) - sum( lgamma(alpha2) ) + sum( (alpha2 - 1) * y_log )
      S_phi <- S_phi + digamma(phi) - sum( mu2 * digamma(alpha2) ) + sum( mu2 * y_log )
      H_phi <- H_phi + trigamma(phi) - sum( mu2^2 * trigamma(alpha2) )

      ## standard full-D softmax Jacobian dmu/deta, then restrict/renormalise
      ## it to the survivor sub-simplex: mu2 = mu[Pi] / sum(mu[Pi])
      J <- matrix(0, nrow = D, ncol = p)
      mu_sub <- mu_i[2:D]
      for ( d in 1:D ) {
        for ( k in 1:p ) {
          J[d, k] <- mu_i[d] * ( (d == (k + 1)) - mu_sub[k] )
        }
      }
      Jp <- J[Pi, , drop = FALSE]                    # |Pi| x p
      Jsum <- Rfast::colsums(Jp)                     # p
      J2 <- ( Jp - outer(mu2, Jsum) ) / Mi            # |Pi| x p, d mu2/d eta

      h_i <- y_log - digamma(alpha2)                  # |Pi|
      psi1_alpha2 <- trigamma(alpha2)
      C <- diag(psi1_alpha2, nrow = length(Pi)) - psi1_phi
      JtCJ <- crossprod(J2, C) %*% J2
      Jt_h <- crossprod(J2, h_i)

      S_i <- as.vector( phi * tcrossprod(x_i, Jt_h) )
      I_i <- phi^2 * kronecker(JtCJ, txi[[ i ]] )

      S_vec <- S_vec + S_i
      I_mat <- I_mat + I_i
    }

    if ( abs(loglik - loglik_old) < tol ) break

    beta_vec <- as.vector(beta)
    beta_vec_new <- beta_vec + solve(I_mat, S_vec)
    beta <- matrix(beta_vec_new, nrow = K, ncol = p)

    phi_new <- phi - S_phi / H_phi
    if ( phi_new <= 0 ) phi_new <- 1e-4
    phi <- phi_new

    loglik_old <- loglik
  }

  colnames(beta) <- paste0("Y", 2:D)
  rownames(beta) <- colnames(x)

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., as.data.frame(xnew) )
    if ( !con )  xnew <- xnew[, -1, drop = FALSE]
    ma <- cbind(1, exp( xnew %*% be ) )
    est <- ma / Rfast::rowsums(ma)  ## fitted values
    colnames(est) <- colnames(y)
  }

  list( iters = iter, loglik = loglik + const, const = const, phi = phi, be = beta, est = est )
}
