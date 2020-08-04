#before running the code: make sure rtools is properly installed on your computer if you use windows OS.
#install.packages("pracma", "statmod", "Matrix", "SparseGrid", "devtools", "R6", "compiler", "randtoolbox", "rngWELL")
#devtools::install_github("TillF/ppso")
library(ppso)

######################################################
##########same elementary designs, fix 0, 12##########
######################################################
rm(list=ls())
# objective function
dcri2<-function(d){
  d<-sort(d)
  #replace to your own path
  setwd("/Users/yushi/Desktop/")
  source("app2/model.R")
  
  # define longitudinal binary model
  model <- Model$new(
    parameter_function = function(mu, b) list(base=mu[1]+b[1], slp=mu[2]+b[2], eff=mu[3]),
    log_likelihood_function = function(y, design, base, slp, eff) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1+eff*trt)*time))), design$time, design$trt)
      log(ifelse(y==1,p,1-p))
    }, 
    simulation_function = function(design, base, slp, eff) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1+eff*trt)*time))), design$time, design$trt)
      as.numeric(p > runif(nrow(design)))
    },
    inverse_simulation_function = function(design, urand, base, slp, eff) {
      if(is.null(urand)) return(seq_along(design$time))
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1+eff*trt)*time))), design$time, design$trt)
      qbinom(urand, 1, prob = p)
    },
    mu = c(-2, 0.09, 5),
    omega = diag(c(0.7^2, 0.17^2)))
  
  # define settings
  settings <- defaults.agq(gq.quad_points = 7,  y_integration.method = "qrmc", y_integration.n_samples = 1000)
  # define design 
  d2<-c(0, d, 12)
  design <- data.frame(time=d2, trt=1)
  # use only half of the samples for each group
  settings$y_integration.n_samples <- ceil(settings$y_integration.n_samples/2)
  # calculate fim for treatment group
  fim_trt <- calc_fim(model, design, settings)
  # calculate fim for placebo group
  fim_plb <- calc_fim(model, transform(design, trt=0), settings)
  # add fims together
  fim <-   50*(fim_trt + fim_plb)
  dcri<-det(fim)^(-1/5) #transform to minimization
  return(dcri)
  rm(list=ls())
}

# initial time
initime2=c(4, 8)
number_of_parameters<-2
number_of_particles <- 40
max_number_of_iterations<-100
max_number_function_calls<-100
parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters))
ptm <- proc.time()
result<-optim_pso(objective_function = dcri2, number_of_parameters = 2, number_of_particles = 40, 
                  max_number_of_iterations = 100,max_number_function_calls= 100, 
                  parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters)), 
                  initial_estimates=initime2, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                  lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm
1/result$value
result


######################################################
##########same elementary designs, flexible###########
######################################################
rm(list=ls())
# objective function
dcri<-function(d){
  d<-sort(d)
  #replace to your own path
  setwd("/Users/yushi/Desktop/")
  source("app2/model.R")
  
  # define longitudinal binary model
  model <- Model$new(
    parameter_function = function(mu, b) list(base=mu[1]+b[1], slp=mu[2]+b[2], eff=mu[3]),
    log_likelihood_function = function(y, design, base, slp, eff) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1+eff*trt)*time))), design$time, design$trt)
      log(ifelse(y==1,p,1-p))
    }, 
    simulation_function = function(design, base, slp, eff) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1+eff*trt)*time))), design$time, design$trt)
      as.numeric(p > runif(nrow(design)))
    },
    inverse_simulation_function = function(design, urand, base, slp, eff) {
      if(is.null(urand)) return(seq_along(design$time))
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1+eff*trt)*time))), design$time, design$trt)
      qbinom(urand, 1, prob = p)
    },
    mu = c(-2, 0.09, 5),
    omega = diag(c(0.7^2, 0.17^2)))
  
  # define settings
  settings <- defaults.agq(gq.quad_points = 7,  y_integration.method = "qrmc", y_integration.n_samples = 1000)
  # define design 
  design <- data.frame(time=d, trt=1)
  # use only half of the samples for each group
  settings$y_integration.n_samples <- ceil(settings$y_integration.n_samples/2)
  # calculate fim for treatment group
  fim_trt <- calc_fim(model, design, settings)
  # calculate fim for placebo group
  fim_plb <- calc_fim(model, transform(design, trt=0), settings)
  # add fims together
  fim <-   50*(fim_trt + fim_plb)
  dcri<-det(fim)^(-1/5) #transform to minimization
  return(dcri)
  rm(list=ls())
}

initime=c(0, 1.5, 6.2, 12)
number_of_parameters<-4
number_of_particles <- 40
max_number_of_iterations<-100
max_number_function_calls<-100
parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters))
ptm <- proc.time()
result<-optim_pso(objective_function = dcri, number_of_parameters = 4, number_of_particles = 40, 
                  max_number_of_iterations = 100,max_number_function_calls= 100,
                  parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters)), 
                  initial_estimates=initime, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                  lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm
1/result$value
result


######################################################
#####different elementary designs, fixed 0, 12########
######################################################
rm(list=ls())
# objective function
dcri.ele2<-function(d){
  d1<-sort(d[1:2])
  d2<-sort(d[3:4])
  #replace to your own path
  setwd("/Users/yushi/Desktop/")
  source("app2/model.R")
  
  # define longitudinal binary model
  model <- Model$new(
    parameter_function = function(mu, b) list(base=mu[1]+b[1], slp=mu[2]+b[2], eff=mu[3]),
    log_likelihood_function = function(y, design, base, slp, eff) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1+eff*trt)*time))), design$time, design$trt)
      log(ifelse(y==1,p,1-p))
    }, 
    simulation_function = function(design, base, slp, eff) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1+eff*trt)*time))), design$time, design$trt)
      as.numeric(p > runif(nrow(design)))
    },
    inverse_simulation_function = function(design, urand, base, slp, eff) {
      if(is.null(urand)) return(seq_along(design$time))
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1+eff*trt)*time))), design$time, design$trt)
      qbinom(urand, 1, prob = p)
    },
    mu = c(-2, 0.09, 5),
    omega = diag(c(0.7^2, 0.17^2)))
  
  # define settings
  settings <- defaults.agq(gq.quad_points = 7,  y_integration.method = "qrmc", y_integration.n_samples = 1000)
  # define design 
  d11<-c(0, d1, 12)
  design.t <- data.frame(time=d11, trt=1)
  d22<-c(0, d2, 12)
  design.p <- data.frame(time=d22, trt=0)
  # use only half of the samples for each group
  settings$y_integration.n_samples <- ceil(settings$y_integration.n_samples/2)
  # calculate fim for treatment group
  fim_trt <- calc_fim(model, design.t, settings)
  # calculate fim for placebo group
  fim_plb <- calc_fim(model, design.p, settings)
  # add fims together
  fim <-   50*(fim_trt + fim_plb)
  dcri<-det(fim)^(-1/5) #transform to minimization
  return(dcri)
  rm(list=ls())
}

initime2=c(1.5, 6.2, 1.5, 6.2)
number_of_parameters<-4
number_of_particles <- 40
max_number_of_iterations<-100
max_number_function_calls<-100
parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters))
ptm <- proc.time()
result<-optim_pso(objective_function = dcri.ele2, number_of_parameters = 4, number_of_particles = 40, 
                  max_number_of_iterations = 100,max_number_function_calls= 100, 
                  parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters)), 
                  initial_estimates=initime2, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                  lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm
1/result$value
result

######################################################
########different elementary designs, flexible########
######################################################
rm(list=ls())
# objective function
dcri.ele<-function(d){
  d1<-sort(d[1:4])
  d2<-sort(d[5:8])
  #replace to your own path
  setwd("/Users/yushi/Desktop/")
  source("app2/model.R")
  
  # define longitudinal binary model
  model <- Model$new(
    parameter_function = function(mu, b) list(base=mu[1]+b[1], slp=mu[2]+b[2], eff=mu[3]),
    log_likelihood_function = function(y, design, base, slp, eff) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1+eff*trt)*time))), design$time, design$trt)
      log(ifelse(y==1,p,1-p))
    }, 
    simulation_function = function(design, base, slp, eff) {
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1+eff*trt)*time))), design$time, design$trt)
      as.numeric(p > runif(nrow(design)))
    },
    inverse_simulation_function = function(design, urand, base, slp, eff) {
      if(is.null(urand)) return(seq_along(design$time))
      p <- mapply(function(time, trt) 1/(1+exp(-(base+slp*(1+eff*trt)*time))), design$time, design$trt)
      qbinom(urand, 1, prob = p)
    },
    mu = c(-2, 0.09, 5),
    omega = diag(c(0.7^2, 0.17^2)))
  
  # define settings
  settings <- defaults.agq(gq.quad_points = 7,  y_integration.method = "qrmc", y_integration.n_samples = 1000)
  # define design 
  design.t <- data.frame(time=d1, trt=1)
  design.p <- data.frame(time=d2, trt=0)
  # use only half of the samples for each group
  settings$y_integration.n_samples <- ceil(settings$y_integration.n_samples/2)
  # calculate fim for treatment group
  fim_trt <- calc_fim(model, design.t, settings)
  # calculate fim for placebo group
  fim_plb <- calc_fim(model, design.p, settings)
  # add fims together
  fim <-   50*(fim_trt + fim_plb)
  dcri<-det(fim)^(-1/5) #transform to minimization
  return(dcri)
  rm(list=ls())
}

initime=c(0, 1.5, 2.8, 12, 0, 4.2, 9.1, 12)
number_of_parameters<-8
number_of_particles <- 40
max_number_of_iterations<-100
max_number_function_calls<-100
parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters))
ptm <- proc.time()
result<-optim_pso(objective_function = dcri.ele, number_of_parameters = 8, number_of_particles = 40, 
                  max_number_of_iterations = 100,max_number_function_calls= 100, 
                  parameter_bounds = cbind(rep(0, number_of_parameters),rep(12, number_of_parameters)), 
                  initial_estimates=initime, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                  lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm
1/result$value
result
