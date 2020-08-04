ParLogLike <- function(y, design, eta, mu = .mu, sqrt.omega=.sqrt.omega.matrix){
  do.call(LogLike, c(list(y=y, design=design),Par(mu, sqrt.omega %*% t(eta))))
}

ParSim <- function(design, eta, mu = .mu, sqrt.omega=.sqrt.omega.matrix){
  do.call(Sim, c(list(design=design),Par(mu, sqrt.omega %*% t(eta))))
}