# Function to estimate Dirichlet precision parameter (alpha_0) using Newton-Raphson
diria0.est <- function(x, tol = 1e-6) {
  n <- dim(x)[1]   ;   D <- dim(x)[2]
  a <- 0   ;   ea <- 1
  slx <- sum( Rfast::Log(x) )
  lik1 <- n * ( lgamma(D * ea) - D * lgamma(ea) ) + (ea - 1) * slx
  grad <- n * D * ea * ( digamma(D * ea) - digamma(ea) ) + ea * slx
  hess <- n * trigamma(D * ea) * (D * ea)^2 - n * D * trigamma(ea) * ea^2 + grad
  a <- a - grad/hess
  ea <- exp(a)
  lik2 <- n * ( lgamma(D * ea) - D * lgamma(ea) ) + (ea - 1) * slx
  i <- 1

  while ( lik2 - lik1 > tol ) {
    i <- i + 1
    lik1 <- lik2
    grad <- n * D * ea * ( digamma(D * ea) - digamma(ea) ) + ea * slx
    hess <- n * trigamma(D * ea) * (D * ea)^2 - n * D * trigamma(ea) * ea^2 + grad
    a <- a - grad/hess
    ea <- exp(a)
    lik2 <- n * ( lgamma(D * ea) - D * lgamma(ea) ) + (ea - 1) * slx
  }
  list(iter = i, loglik = lik2, a0 = ea)
}
