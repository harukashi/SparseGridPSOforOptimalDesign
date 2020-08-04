library(pracma)
library(statmod)
library(Matrix)
library(SparseGrid)

integrate_gq <- function(fun, dimensions, center=rep(0,dimensions),
                         scale=diag(1,dimensions), settings=filter_settings(defaults.agq(), "gq"), ...){
  additional.args <- c(list(), list(...))
  n.quad.points <- settings$quad_points
  if(is.character(dimensions)) {
    ndim <- length(dimensions)
  }else{
    ndim <- dimensions
  }
  if(dimensions==0) {
    n.quad.points <- 1
    ndim <- 1
  }

  sg <- createSparseGrid("KPN", ndim, k = 3)
  sqrt_scale <- t(chol(nearPD(scale)$mat))
  adp_grid_points <- matrix(apply(sg$nodes, 1, function(z) as.matrix(center+sqrt_scale%*%z)), ncol = ndim, byrow = T)
  adp_weights <- sg$weights * apply(adp_grid_points, 1, function(z) prod(dnorm(z)))/apply(sg$nodes, 1, function(z)prod(dnorm(z)))
  call.args <- c(additional.args, list(adp_grid_points))
  grid.results <- do.call(fun, call.args)
  det(sqrt_scale) * grid.results %*% adp_weights 
}




integrate_mc <- function(fun, dimensions, settings=filter_settings(defaults.agq(), "mc"), ...){
  additional.args <- c(list(), list(...))
  n.samples <- settings$n_samples
  if(is.character(dimensions)) {
    ndim <- length(dimensions)
  }else{
    ndim <- dimensions
  }
  if(ndim>0){ 
    param.samples <- matrix(rnorm(n.samples*ndim, 0, 1), ncol=ndim) 
  }else{
    param.samples <- matrix(nrow=n.samples)
  }
  if(is.character(dimensions)) {
    colnames(param.samples) <- dimensions
  }
  call.args <- c(additional.args, list(param.samples))  
  mc.results <- do.call(fun, call.args)
  if(!is.matrix(mc.results)) return(mean(mc.results))
  return(rowMeans(mc.results))
}

integrate_mc_lhs <- function(fun, dimensions, settings=filter_settings(defaults.agq(), "mc"), ...){
  require(lhs)
  additional.args <- c(list(), list(...))
  n.samples <- settings$n_samples
  if(is.character(dimensions)) {
    ndim <- length(dimensions)
  }else{
    ndim <- dimensions
  }
  
  inv_samples <- randomLHS(n.samples,dimensions)
  call.args <- c(additional.args, list(inv_samples))  
  mc.results <- do.call(fun, call.args)
  if(!is.matrix(mc.results)) return(mean(mc.results))
  return(rowMeans(mc.results))
}

get_sobol_closure <- function(){
  .dim <- 0
  function(samples, dim){
    if(.dim==dim) sobol(samples, dim = dim, scrambling = T, init=F)
    else{
      .dim <<- dim
      sobol(samples, dim = dim, scrambling = T, init=T)
    }
  }
}

my_sobol <- get_sobol_closure()

integrate_qrmc <- function(fun, dimensions, settings=filter_settings(defaults.agq(), "qrmc"), ...){
  require(randtoolbox)
  additional.args <- c(list(), list(...))
  n.samples <- settings$n_samples
  if(is.character(dimensions)) {
    ndim <- length(dimensions)
  }else{
    ndim <- dimensions
  }

  inv_samples <- my_sobol(n.samples, dimensions)
  if(dimensions==1) inv_samples <- matrix(inv_samples, ncol=1)
  call.args <- c(additional.args, list(inv_samples))  
  mc.results <- do.call(fun, call.args)
  if(!is.matrix(mc.results)) return(mean(mc.results))
  return(rowMeans(mc.results))
}
