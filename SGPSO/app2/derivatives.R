
h <- 1e-6
JacobianPar <- function(mu = .mu, b=rep(0, .n.b)){
  # Calculates the Jacbian of the parameter function wrt mu and b
  # 
  # Args:
  #  mu: parameter vector of fixed effects to take the derivative at
  #  b:  parameter vector of random effects to take the derivative at
  #
  # Returns:
  #  Jacobian (n.parX(n.mu+n.b)) wrt mu and to b 
  par0 <- unlist(Par(mu, b))
  grad.mu <- sapply(seq_along(mu), function(i) (unlist(Par(mu+h*(seq_along(mu)==i), b))-par0))/h
  grad.b <- sapply(seq_along(b), function(i) (unlist(Par(mu, b+h*(seq_along(b)==i)))-par0))/h
  cbind(grad.mu, grad.b)
}

GradientLogLike <- function(y, design, ...){
  # Calculates the gradient of the log likelihood function wrt to its parameters
  # 
  # Args:
  #  y:       data
  #  design:  design
  #  ...:     parameters to the log-likelihood function to calculate derivative for
  #
  # Returns:
  #  Vector (n.par) of derivatives
  
  params <- list(...)
  ll0 <- do.call(LogLike, c(list(y=y, design=design), params))
  grad <- sapply(seq_along(params), function(i) {
    params[[i]] <- params[[i]]+h
    (do.call(LogLike, c(list(y=y, design=design), params))-ll0)/h
    })
  grad
}

JacobianCholeskyFactor <- function(eta=rep(0, .n.b), omega=.omega.matrix){
  # Calculates the Jacobian of the Cholseky factor wrt to the diagonals of omega
  # 
  # Args:
  #  eta:       vector (n.b) of standardized random effects
  #  omega:     omega matrix
  #
  # Returns:
  #  Matrix (n.bXn.b) with derivatives of elements of b (rows) wrt to the diagnoal elements of omega (columns)
  
  # Use nice formula by Simo et al. for partial derivative of Cholesky factorization (works with general omega)
  # TODO: move outside of function to speed up calculation (omega is constant)
  L <- t(chol(omega))
  L.inv <- solve(L)
  domega <- omega
  domega[] <- 0
  dL <- sapply(seq_len(nrow(omega)), function(i) {
    domega[i,i] <- 1
    M <- L%*%domega%*%t(L.inv)
    M[upper.tri(M)] <- 0
    diag(M) <- 0.5*diag(M)
    L%*%M
  })
  dL
}

tr <- function(m) sum(diag(m))
get_jacobian_pardensity_wrt_omega <- function(eta){
  inv_sqrt_omega <- solve(.sqrt.omega.matrix)
  domega <- .omega.matrix
  domega[] <- 0
  d_sqrt_omega_d_omega <- sapply(seq_len(nrow(.omega.matrix)), function(i) {
    domega[i,i] <- 1
    M <- inv_sqrt_omega%*%domega%*%t(inv_sqrt_omega)
    M[upper.tri(M)] <- 0
    diag(M) <- 0.5*diag(M)
    .sqrt.omega.matrix%*%M
  })
  
  (-det(inv_sqrt_omega)*tr(d_sqrt_omega_d_omega%*%inv_sqrt_omega)+det(inv_sqrt_omega)%*%eta%*%inv_sqrt_omega%*%d_sqrt_omega_d_omega%*%eta)
}

GradientComposite <- function(y, design, eta, mu=.mu, sqrt.omega=.sqrt.omega.matrix){
  # Calculates the gradient of the composite log likelihood function wrt to mu and b
  # 
  # Args:
  #  y:       data
  #  design:  design
  #  b:  parameter vector of random effects to take the derivative at
  #  mu: parameter vector of fixed effects to take the derivative at
  #
  # Returns:
  #  Vector (n.mu+n.b) of derivatives
  inv_sqrt_omega <- solve(.sqrt.omega.matrix)
  b <- sqrt.omega %*% t(eta)
  dLLdP <- do.call(GradientLogLike, c(list(y, design), Par(mu, b)))
  dPdMuB <- JacobianPar(mu, b)
  dLLdMuB <- dLLdP %*% dPdMuB
  dLLdMu <- dLLdMuB[seq_len(.n.mu)] 
  dLLdOmega <- get_jacobian_pardensity_wrt_omega(eta)
  like <- exp(ParLogLike(y, design, eta, mu))
  return(c(dLLdMu*det(inv_sqrt_omega), dLLdOmega)*like)
}

