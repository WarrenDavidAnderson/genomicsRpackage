

#' System identification
#'
#' This function performs system identification to infer the connectivity and dynamics of a 
#' linear network model (Anderson et al., 2017).
#' @param datHMF The data matrix with scaled measurements measurements. Scaled center measurements are provided for 
#' each measurement/annotation (e.g., gene/organ) combination. This format is produced by scale_zeroOne(). Defaults to datHMF.
#' @param N_mRanges The number of m-value sets to evaluate.
#' @param alpha The L1 regularization terms for constraining the regression. Defaults to seq(0,1,0.2).
#' @param nlambda The number of L2 regularization terms for consideration. Defaults to 10.
#' @param Mall A set of m-values. Defaults to NULL.
#' @param Msc A scale factor that determines the M value (M = max(m)) corresponding to the largest range of m-values. 
#' Mmax <- ceiling(Msc * Np / 2) where Np = number of parameters.
#' @param n Highest order of the system (n=1 for the first degree system that this analysis was designed for).
#' @param dt Simulation time step
#' @param epsilon Term for computing simulation error. Defaults to 10^(-10).
#' @param Harg0 If omega_0 is equal to zero (x=0), set omega_0 to Harg0 (w=Harg0). 
#' @return Error data for all simulations, simulation data for the best simulation, a parameter matrix for the best simulation.
#' @export
#' @examples err0 = HMF_fit(datHMF,N_mRanges=2,dt=0.1,alpha=c(0.2,0.6))
#' @references \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005627}
HMF_fit = function(datHMF=NULL,N_mRanges=10,alpha=NULL,nlambda=NULL,Mall=NULL,Msc=2,n=1,dt=NULL,epsilon=10^(-10),Harg0=0.000001){
  
  # regularization terms
  if(is.null(alpha)==F){alph = alpha}
  if(is.null(alpha)==T){alph  = seq(0,1,0.2)}
  if(is.null(nlambda)==F){nlambda = nlambda}
  if(is.null(nlambda)==T){nlambda  = 10}
  
  # total number of model parameters
  numPar <- ((ncol(datHMF)-2)^2)*(length(unique(datHMF[,2]))^2) 
  
  # model initial conditions for simulation
  ic = get_init_condits(datHMF)
  
  # X and Y data, dims: gene, organ, m value (in Mall)
  mvals = get_M_vals(datHMF=datHMF,Msc=Msc)
  Mall = mvals$Mall
  XY = getXY(datHMF,Mall=Mall,Msc=Msc,n=n,Harg0=Harg0)
  
  mspecs = get_M_indices(N_mRanges) # m-value indices
  Mspectrum = mspecs$Mspectrum # indices in the vector of m-values 
  indM = mspecs$indM # select which m-value ranges to plot
  
  # loop through alpha
  simulation_error = list()
  parameters = list()
  for (alp in 1:length(alph)) {
    
    # space for saving rate constant data
    fitCoefs <- array(c(0),c(numPar,nlambda,length(indM))) # fit coefficients
    MLAM <- matrix(c(0),nlambda*length(indM),5)
    
    # loop through every set of m-value ranges
    # implement parameter estimation--regression for each set of m-values
    for (mm in 1:length(indM)) {
      
      # m-indices for a given range of m-values
      mRange = get_M_indices(N_mRanges,mvals,Msc,indM[mm])$mRange
      
      # get matrices for regression
      Yreg = get_Yreg(datHMF,XY$Y,numPar,mRange)
      Xreg = get_Xreg(datHMF,XY$X,numPar,mRange)
      
      # fit the regression model using regularization
      # alpha <- alph[alp]
      fit <- glmnet::glmnet( data.matrix(Xreg), data.matrix(Yreg), family="gaussian",
                     alpha=alph[alp], intercept=F, standardize=F, nlambda=nlambda)
      rateConst = get_rateConst(fit$beta,nlambda)
      #fitCoefs[,1:ncol(rateConst),mm] <- rateConst
      indLam <- ((mm-1)*nlambda+1) : (mm*nlambda)
      MLAM[indLam,1] <- mm # m-value index (1,2,3,...)
      MLAM[indLam,2] <- c(1:nlambda) # lambda index (1,2,3,...)
      MLAM[indLam,3] <- indM[mm] # which value from indM
      MLAM[indLam,4] <- fit$lambda # numeric lambda value
      colnames(MLAM) = c("m-value_index","lambda_index","imdM_value","lambda_value","error")
      
      # model simulation - get errors
      tlen = max(unique(datHMF[,1])) - min(unique(datHMF[,1]))
      Nr = length(unique(datHMF[,2])) # number of organs (annotations)
      MLAM[indLam,5] <- sim_error(datHMF,Nr,tlen,rateConst,ic,dt,epsilon,numPar)
      
    } # end mm loop (m-value ranges)
    
    # paraOut = list()
    # for (mv in 1:length(indM)) {
    #   paraOut[[mv]] <- fitCoefs[,,mv]
    # } 
    # parameters[[alp]] = paraOut
    simulation_error[[alp]] = MLAM
    
  } # end alp loop
  
  # select best simulation 
  simulate_error = format_error(simulation_error,alph)
  best_params = select_best(simulate_error)
  
  # simulate the best model
  # best_params: col1=alpha, col2=mm, col4=indM[mm], col5=lambda
  indM_best = best_params$imdM_value
  mRange_best = get_M_indices(N_mRanges,mvals,Msc,indM_best)$mRange
  tsim <- seq(0,tlen,by=dt)
  best_sims = sim_best(best_params,mRange_best,tsim,numPar,XY,datHMF,ic)
  best_sims$simulations[,1] = best_sims$simulations[,1] + min(datHMF[,1])
  
  return(list(errors = simulate_error, best_sim = best_sims$simulations, best_params = best_sims$param_matrix))
  
} # end HMF_fit()