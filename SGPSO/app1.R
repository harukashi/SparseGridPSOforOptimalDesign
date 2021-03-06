#install.packages("devtools")
#devtools::install_github("TillF/ppso")
library(ppso)

funFIMem <- function(equation,paramName,beta,o,sigma,t_group,Trand,d,PropSubjects,nbTot){
  #Name of the fixed effect parameters and sampling times
  paramF<-c(paramName,"t")
  
  #model equation
  form1<- equation
  
  PSI<-c(paramName,"sig.inter","sig.slope")
  lpsi<-length(PSI)
  
  
  #(Diagonal Matrix of) variance for inter-subject random effects:
  omega<-diag(o)
  
  
  #Random effect model (1) = additive  (2) = exponential 
  #------------------------------------------------------------------
  
  if ( Trand == 1 ) {
    form11 <- form1
    for(i in 1:length(paramName)){
      form11 <- gsub(paramName[i], paste0("(",paramName[i],"+b)"), form11)
      
    }
    
  } else {
    form11 <- form1
    for(i in 1:length(paramName)){
      form11 <- gsub(paramName[i], paste0("(",paramName[i],"*exp(b))"), form11)
    }
  }
  
  #gather all groups of protocol
  M_f<-list()
  M_F <- matrix(rep(0),nrow=lpsi+length(paramName),ncol=lpsi+length(paramName))
  for(q in 1:length(t_group)){
    t<-c()
    t<-c(t_group[[q]])
    
    #dose value
    dose<-d[q]
    
    #calculate matrix E for n individuals
    equatf <- parse(text = form11, n=-1)
    f<-function(paramF){eval(equatf[[1]])}
    
    #Fixed effects parameters values
    for(i in 1:length(paramName)){
      assign(paramName[i],beta[i])
    }
    
    param <- c(beta,t)
    #calculate the predictions with fixed effects
    b <- 0
    fixed<-f(param)
    
    #Standard deviation of residual error (sig.inter+sig.slope*f)^2:
    sig.inter<-sigma[1]
    sig.slope<-sigma[2]
    
    var<-diag((sig.inter+sig.slope*fixed)^2)
    
    
    #calculate variance Var
    form2 <- paste0("( sig.inter + sig.slope * fixed )^2")
    Vmodel<- parse(text = form2)
    
    #get derivatives for fixed parameters
    df<-deriv(equatf[[1]],PSI)
    mdf<-attributes(eval(df))$gradient
    #delete the last two columns (correspond to sig.inter and sig.slope) 
    mdfi <- mdf[,-c(length(PSI)-1,length(PSI))]
    #complete derivative for exponential random effect model 
    if(Trand ==2 ){
      mdfie <- mdfi %*% diag(beta)
    }else {mdfie <- mdfi}
    
    
    #calculate variance Vi
    Vi <- mdfie %*% omega %*% t(mdfie) + var
    
    #get derivatives of sigma
    dv<-deriv(Vmodel[[1]],PSI)
    mdv<-attributes(eval(dv))$gradient
    
    
    #calculate matrix part A
    M_A <- t(mdfi) %*% solve(Vi) %*% mdfi
    
    #complete the rest of the matrix with 0
    for(i in 1:length(PSI)){
      M_A <- cbind(M_A,0)
      M_A <- rbind(M_A,0)
    }
    
    #calculate matrix part B
    #initialize the matrix with 0
    M_B <- matrix(rep(0),nrow=length(PSI)+length(paramName),ncol=length(PSI)+length(paramName))
    #prepare a list of matrix of derivatives of sigma to simplify usage
    m<-list()
    for(i in (length(paramName)+1):lpsi){
      if(length(t)==1){
        m[[i]] <- mdv[i]
      }else{
        m[[i]] <- diag(mdv[,i])
      }
      
    }
    #calculate first three rows of part B
    for(i in 1:length(paramName)){
      
      for(j in 1:length(paramName)){
        M_B[length(paramName)+i,length(paramName)+j] <- 1/2 * sum(diag(((mdfie[,i] %*% t(mdfie[,i])) %*% solve(Vi) %*% (mdfie[,j] %*% t(mdfie[,j])) %*% solve(Vi))))
      }
      for(j in (length(PSI)-1):length(PSI)){
        M_B[length(paramName)+i,length(paramName)+j] <- 1/2 * sum(diag(((mdfie[,i] %*% t(mdfie[,i])) %*% solve(Vi) %*% m[[j]] %*% solve(Vi))))
      }
      
    }
    #calculate the last two rows of partB
    for(i in (length(PSI)-1):length(PSI)){
      for(j in 1:length(paramName)){
        M_B[length(paramName)+i,length(paramName)+j] <- 1/2 * sum(diag( m[[i]] %*% solve(Vi) %*% (mdfie[,j] %*% t(mdfie[,j])) %*% solve(Vi)))
      }
      for(j in (length(PSI)-1):length(PSI)){
        M_B[length(paramName)+i,length(paramName)+j] <- 1/2 * sum(diag(m[[i]] %*% solve(Vi) %*% m[[j]] %*% solve(Vi)))
      }
    }
    
    M_f[[q]] <- (M_A+M_B)*PropSubjects[q]
    M_F <-M_F+M_f[[q]]
  }
  M_F <- M_F *nbTot
  
  #set names for vectors 
  fname<-c()
  for(n in 1:length(paramName)){
    fname<-c(fname,paste0("u_",paramName[n]))
  }
  for(n in 1:length(paramName)){
    fname<-c(fname,paste0("w_",paramName[n]))
  }
  fname<-c(fname,"sig.inter","sig.slope")
  rownames(M_F) <- fname
  colnames(M_F) <- fname
  
  if(sig.slope ==0){
    M_F <- M_F[,-c(lpsi+length(paramName))]
    M_F <- M_F[-c(lpsi+length(paramName)),]
    PSI <- PSI[-c(lpsi)]
    M_F<-M_F[,-10]
    M_F<-M_F[-10,]
  }
  if(sig.inter == 0){
    M_F <- M_F[,-c(lpsi+length(paramName)-1)]
    M_F <- M_F[-c(lpsi+length(paramName)-1),]
    PSI <- PSI[-c(lpsi-1)]
  }
  if(length(t)==1){
    return(M_F)
  }else{
    deterFim <- det(M_F)
    if(is.na(deterFim)==FALSE&deterFim>0)
    {CritereDopt <- deterFim^(1/10)}
    else{CritereDopt<-0}
    
    return(CritereDopt)
    
  }
}


#####################################################
#################Fix sampling times##################
#####################################################
Dopt1<-function(time){
  res<-funFIMem("exp(VA0) + (1-exp(-exp(Kpr)*t))* (Emax*dose/(exp(ED50)+dose) - exp(Beta)*exp(VA0))",c("VA0","Kpr","Beta", "Emax", "ED50"),
                c(log(55),log(0.005),log(0.2),30,log(150)),c(0.07,0.5,1,150,0),c(sqrt(28), 0),
                list(c(0, 224, 420, 672), c(0, 224, 420, 672), c(0, 224, 420, 672), c(0, 224, 420, 672)),1,c(time[1], time[2], time[3], time[4]),
                c(time[5]/sum(time[5:8]), time[6]/sum(time[5:8]), time[7]/sum(time[5:8]), time[8]/sum(time[5:8])),300)
  return(-res)
}

ptm <- proc.time()
number_of_parameters<-8
parameter_bounds = cbind(rep(0, number_of_parameters),rep(500, number_of_parameters))
result1<-optim_pso(objective_function = Dopt1, number_of_parameters = 8, number_of_particles = 40, 
                   max_number_of_iterations = 1000, max_number_function_calls= 1000,
                   parameter_bounds = cbind(rep(0, number_of_parameters),rep(500, number_of_parameters)), 
                   initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                   lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

doses<-result1$par[1:4]
allocation1<-round(result1$par[5]/sum(result1$par[5:8])*300, 0)
allocation2<-round(result1$par[6]/sum(result1$par[5:8])*300, 0)
allocation3<-round(result1$par[7]/sum(result1$par[5:8])*300, 0)
allocation4<-round(result1$par[8]/sum(result1$par[5:8])*300, 0)

result1$value
doses
c(allocation1, allocation2, allocation3, allocation4)

#####################################################
####################Fix allocation###################
#####################################################
Dopt2<-function(time){
  res<-funFIMem("exp(VA0) + (1-exp(-exp(Kpr)*t))* (Emax*dose/(exp(ED50)+dose) - exp(Beta)*exp(VA0))",c("VA0","Kpr","Beta", "Emax", "ED50"),
                c(log(55),log(0.005),log(0.2),30,log(150)),c(0.07,0.5,1,150,0),c(sqrt(28), 0),
                list(time[1:4],time[5:8],time[9:12],time[13:16]),1,c(time[17], time[18], time[19], time[20]),c(0.25, 0.25, 0.25, 0.25),300)
  return(-res)
}

ptm <- proc.time()
number_of_parameters<-20
parameter_bounds = cbind(rep(0, number_of_parameters),rep(500, number_of_parameters))
result2<-optim_pso(objective_function = Dopt2, number_of_parameters = 20, number_of_particles = 40, 
                   max_number_of_iterations = 1000,max_number_function_calls= 1000,
                   parameter_bounds = cbind(rep(0, number_of_parameters),c(rep(672, 16), rep(500, 4))), 
                   initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                   lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-sort(round(result2$par[1:4],2))
t_group2<-sort(round(result2$par[5:8],2))
t_group3<-sort(round(result2$par[9:12],2))
t_group4<-sort(round(result2$par[13:16],2))
doses<-round(result2$par[17:20], 2)

result2$value
doses
t_group1
t_group2
t_group3
t_group4


#####################################################
#######################Fix doses#####################
#####################################################
Dopt3<-function(time){
  res<-funFIMem("exp(VA0) + (1-exp(-exp(Kpr)*t))* (Emax*dose/(exp(ED50)+dose) - exp(Beta)*exp(VA0))",c("VA0","Kpr","Beta", "Emax", "ED50"),
                c(log(55),log(0.005),log(0.2),30,log(150)),c(0.07,0.5,1,150,0),c(sqrt(28), 0),
                list(time[1:4],time[5:8],time[9:12],time[13:16]),1,c(0, 150, 300, 500),c(time[17]/sum(time[17:20]), time[18]/sum(time[17:20]), time[19]/sum(time[17:20]), time[20]/sum(time[17:20])),300)
  return(-res)
}

number_of_parameters<-20
parameter_bounds = cbind(rep(0, number_of_parameters),rep(672, number_of_parameters))
ptm <- proc.time()
result3<-optim_pso(objective_function = Dopt3, number_of_parameters = 20, number_of_particles = 40, 
                   max_number_of_iterations = 1000,max_number_function_calls= 1000, 
                   parameter_bounds = cbind(rep(0, number_of_parameters),rep(672, number_of_parameters)), 
                   initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                   lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-round(result3$par[1:4],2)
t_group2<-round(result3$par[5:8],2)
t_group3<-round(result3$par[9:12],2)
t_group4<-round(result3$par[13:16],2)
allocation1<-round(result3$par[17]/sum(result3$par[17:20])*300, 0)
allocation2<-round(result3$par[18]/sum(result3$par[17:20])*300, 0)
allocation3<-round(result3$par[19]/sum(result3$par[17:20])*300, 0)
allocation4<-round(result3$par[20]/sum(result3$par[17:20])*300, 0)
design.matrix<-matrix(c(t_group1, allocation1, t_group2, allocation2, t_group3, allocation3, t_group4, allocation4), byrow=TRUE, ncol=5)
colnames(design.matrix)<-c("time1", "time2", "time3", "time4", "N patients")

result3$value
design.matrix

#####################################################
##################All flexible#######################
#####################################################
Dopt4<-function(time){
  res<-funFIMem("exp(VA0) + (1-exp(-exp(Kpr)*t))* (Emax*dose/(exp(ED50)+dose) - exp(Beta)*exp(VA0))",c("VA0","Kpr","Beta", "Emax", "ED50"),
                c(log(55),log(0.005),log(0.2),30,log(150)),c(0.07,0.5,1,150,0),c(sqrt(28), 0),
                list(time[1:4],time[5:8],time[9:12],time[13:16]),1,c(time[21], time[22], time[23], time[24]),c(time[17]/sum(time[17:20]), time[18]/sum(time[17:20]), time[19]/sum(time[17:20]), time[20]/sum(time[17:20])),300)
  return(-res)
}

number_of_parameters<-24
parameter_bounds = cbind(rep(0, number_of_parameters),c(rep(672, 20), rep(500, 4)))
ptm <- proc.time()
result4<-optim_pso(objective_function = Dopt4, number_of_parameters = 24, number_of_particles = 40, 
                   max_number_of_iterations = 1000,max_number_function_calls= 1000, 
                   parameter_bounds = cbind(rep(0, number_of_parameters),c(rep(672, 20), rep(500, 4))), 
                   initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                   lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-round(result4$par[1:4],2)
t_group2<-round(result4$par[5:8],2)
t_group3<-round(result4$par[9:12],2)
t_group4<-round(result4$par[13:16],2)
allocation1<-round(result4$par[17]/sum(result4$par[17:20])*300, 0)
allocation2<-round(result4$par[18]/sum(result4$par[17:20])*300, 0)
allocation3<-round(result4$par[19]/sum(result4$par[17:20])*300, 0)
allocation4<-round(result4$par[20]/sum(result4$par[17:20])*300, 0)
dose1<-round(result4$par[21], 1)
dose2<-round(result4$par[22], 1)
dose3<-round(result4$par[23], 1)
dose4<-round(result4$par[24], 1)
design.matrix<-matrix(c(t_group1, allocation1,dose1, t_group2, allocation2,dose2, t_group3, allocation3,dose3, t_group4, allocation4, dose4), byrow=TRUE, ncol=6)
colnames(design.matrix)<-c("time1", "time2", "time3", "time4", "N patients", "dose")

result4$value
design.matrix

