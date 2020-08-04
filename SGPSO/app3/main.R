#before running the code: make sure rtools is properly installed on your computer if you use windows OS.
#install.packages("pracma", "statmod", "Matrix", "SparseGrid", "devtools", "R6", "compiler", "randtoolbox", "rngWELL")
#devtools::install_github("TillF/ppso")
library(ppso)

#############################################
##########all flexible#######################
rm(list=ls())
# objective function
dcri.ele<-function(d){
  d1<-sort(d[1:4])
  d2<-sort(d[5:8])
  #replace to your own path
  setwd("/Users/yushi/Desktop/")
  source("app3/model.R")
  
  model <- Model$new(
    parameter_function = function(mu, b) list(base=mu[1]+b[1], slp=mu[2]+b[2], eff=mu[3]+b[3], theta3 = mu[4]+b[4]),
    log_likelihood_function = function(y, design, base, slp, eff, theta3) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1-eff*trt)*exp(theta3*time-1)))), design$time, design$trt)
      log(ifelse(y==1,p,1-p))
    }, 
    simulation_function = function(design, base, slp, eff, theta3) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1-eff*trt)*exp(theta3*time-1)))), design$time, design$trt)
      as.numeric(p > runif(nrow(design)))
    },
    inverse_simulation_function = function(design, urand, base, slp, eff, theta3) {
      if(is.null(urand)) return(seq_along(design$time))
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1-eff*trt)*exp(theta3*time-1)))), design$time, design$trt)
      qbinom(urand, 1, prob = p)
    },
    mu = c(-1, 4, 0.4, 0.33),
    omega = diag(c(.5, 4, .01, .01)))
  
  settings <- defaults.agq(y_integration.method = "qrmc", y_integration.n_samples = 500)
  design.t <- data.frame(time=d1, trt=1)
  design.p <- data.frame(time=d2, trt=0)
  settings$y_integration.n_samples <- ceil(settings$y_integration.n_samples/2)
  fim_trt <- calc_fim(model, design.t, settings)
  fim_plb <- calc_fim(model, design.p, settings)
  fim <-   50*(fim_trt + fim_plb)
  dcri<-det(fim)^(-1/8) #transform to minimization
  return(dcri)
  rm(list=ls())
}

initime=c(0, 1.4, 2.0, 7.6, 0, 1.4, 2.0, 7.6)
number_of_parameters<-8
number_of_particles <- 40
max_number_of_iterations<-50
max_number_function_calls<-50
parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters))
ptm <- proc.time()
result<-optim_pso(objective_function = dcri.ele, number_of_parameters = 8, number_of_particles = 40, 
                  max_number_of_iterations = 50,max_number_function_calls= 50, 
                  parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters)), 
                  initial_estimates=initime, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                  lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm
1/result$value
result



######################################################
#####fixed 0, 12#####################################
rm(list=ls())
dcri.ele2<-function(d){
  d1<-sort(d[1:2])
  d2<-sort(d[3:4])
  #replace to your own path
  setwd("/Users/yushi/Desktop/")
  source("app3/model.R")
  
  model <- Model$new(
    parameter_function = function(mu, b) list(base=mu[1]+b[1], slp=mu[2]+b[2], eff=mu[3]+b[3], theta3 = mu[4]+b[4]),
    log_likelihood_function = function(y, design, base, slp, eff, theta3) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1-eff*trt)*exp(theta3*time-1)))), design$time, design$trt)
      log(ifelse(y==1,p,1-p))
    }, 
    simulation_function = function(design, base, slp, eff, theta3) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1-eff*trt)*exp(theta3*time-1)))), design$time, design$trt)
      as.numeric(p > runif(nrow(design)))
    },
    inverse_simulation_function = function(design, urand, base, slp, eff, theta3) {
      if(is.null(urand)) return(seq_along(design$time))
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1-eff*trt)*exp(theta3*time-1)))), design$time, design$trt)
      qbinom(urand, 1, prob = p)
    },
    mu = c(-1, 4, 0.4, 0.33),
    omega = diag(c(.5, 4, .01, .01)))
  
  settings <- defaults.agq(y_integration.method = "qrmc", y_integration.n_samples = 500)
  d11<-c(0, d1, 12)
  design.t <- data.frame(time=d11, trt=1)
  d22<-c(0, d2, 12)
  design.p <- data.frame(time=d22, trt=0)
  settings$y_integration.n_samples <- ceil(settings$y_integration.n_samples/2)
  fim_trt <- calc_fim(model, design.t, settings)
  fim_plb <- calc_fim(model, design.p, settings)
  fim <-   50*(fim_trt + fim_plb)
  dcri<-det(fim)^(-1/8) #transform to minimization
  return(dcri)
  rm(list=ls())
}


initime2=c(1.5, 2.8, 4.2, 9.1)
number_of_parameters<-4
number_of_particles <- 40
max_number_of_iterations<-50
max_number_function_calls<-50
parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters))
ptm <- proc.time()
result<-optim_pso(objective_function = dcri.ele2, number_of_parameters = 4, number_of_particles = 40, 
                  max_number_of_iterations = 50,max_number_function_calls= 50, 
                  parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters)), 
                  initial_estimates=initime2, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                  lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm
1/result$value
result



#############################################
##########4 time points (same elementary)####
rm(list=ls())
dcri.ele<-function(d){
  d1<-sort(d)
  #replace to your own path
  setwd("/Users/yushi/Desktop/")
  source("app3/model.R")
  
  model <- Model$new(
    parameter_function = function(mu, b) list(base=mu[1]+b[1], slp=mu[2]+b[2], eff=mu[3]+b[3], theta3 = mu[4]+b[4]),
    log_likelihood_function = function(y, design, base, slp, eff, theta3) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1-eff*trt)*exp(theta3*time-1)))), design$time, design$trt)
      log(ifelse(y==1,p,1-p))
    }, 
    simulation_function = function(design, base, slp, eff, theta3) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1-eff*trt)*exp(theta3*time-1)))), design$time, design$trt)
      as.numeric(p > runif(nrow(design)))
    },
    inverse_simulation_function = function(design, urand, base, slp, eff, theta3) {
      if(is.null(urand)) return(seq_along(design$time))
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1-eff*trt)*exp(theta3*time-1)))), design$time, design$trt)
      qbinom(urand, 1, prob = p)
    },
    mu = c(-1, 4, 0.4, 0.33),
    omega = diag(c(.5, 4, .01, .01)))
  
  settings <- defaults.agq(y_integration.method = "qrmc", y_integration.n_samples = 500)
  design.t <- data.frame(time=d1, trt=1)
  design.p <- data.frame(time=d1, trt=0)
  settings$y_integration.n_samples <- ceil(settings$y_integration.n_samples/2)
  fim_trt <- calc_fim(model, design.t, settings)
  fim_plb <- calc_fim(model, design.p, settings)
  fim <-   50*(fim_trt + fim_plb)
  dcri<-det(fim)^(-1/8) #transform to minimization
  return(dcri)
  rm(list=ls())
}

initime=c(0, 1.4, 2.0, 7.6)
number_of_parameters<-4
number_of_particles <- 40
max_number_of_iterations<-50
max_number_function_calls<-50
parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters))
ptm <- proc.time()
result<-optim_pso(objective_function = dcri.ele, number_of_parameters = 4, number_of_particles = 40, 
                  max_number_of_iterations = 50,max_number_function_calls= 50, 
                  parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters)), 
                  initial_estimates=initime, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                  lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm
1/result$value
result


######################################################
#####fixed 0, 12 (same elementary)####################
rm(list=ls())
dcri.ele2<-function(d){
  d1<-sort(d)
  #replace to your own path
  setwd("/Users/yushi/Desktop/")
  source("app3/model.R")
  
  model <- Model$new(
    parameter_function = function(mu, b) list(base=mu[1]+b[1], slp=mu[2]+b[2], eff=mu[3]+b[3], theta3 = mu[4]+b[4]),
    log_likelihood_function = function(y, design, base, slp, eff, theta3) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1-eff*trt)*exp(theta3*time-1)))), design$time, design$trt)
      log(ifelse(y==1,p,1-p))
    }, 
    simulation_function = function(design, base, slp, eff, theta3) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1-eff*trt)*exp(theta3*time-1)))), design$time, design$trt)
      as.numeric(p > runif(nrow(design)))
    },
    inverse_simulation_function = function(design, urand, base, slp, eff, theta3) {
      if(is.null(urand)) return(seq_along(design$time))
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1-eff*trt)*exp(theta3*time-1)))), design$time, design$trt)
      qbinom(urand, 1, prob = p)
    },
    mu = c(-1, 4, 0.4, 0.33),
    omega = diag(c(.5, 4, .01, .01)))
  
  settings <- defaults.agq(y_integration.method = "qrmc", y_integration.n_samples = 500)
  d11<-c(0, d1, 12)
  design.t <- data.frame(time=d11, trt=1)
  design.p <- data.frame(time=d11, trt=0)
  settings$y_integration.n_samples <- ceil(settings$y_integration.n_samples/2)
  fim_trt <- calc_fim(model, design.t, settings)
  fim_plb <- calc_fim(model, design.p, settings)
  fim <-   50*(fim_trt + fim_plb)
  dcri<-det(fim)^(-1/8) #transform to minimization
  return(dcri)
  rm(list=ls())
}


initime2=c(1.5, 6.2)
number_of_parameters<-2
number_of_particles <- 40
max_number_of_iterations<-50
max_number_function_calls<-50
parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters))
ptm <- proc.time()
result<-optim_pso(objective_function = dcri.ele2, number_of_parameters = 2, number_of_particles = 40, 
                  max_number_of_iterations = 50,max_number_function_calls= 50, 
                  parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters)), 
                  initial_estimates=initime2, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                  lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm
1/result$value
result


