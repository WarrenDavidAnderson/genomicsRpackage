

#' The continuous Hartley transform
#'
#' This function performs linear interpolation between data points and implements the continuous Hartley transform
#' @param w This is omega_0 (Equation 4, Anderson et al., 2017). Defaults to NULL.
#' @param t Vector of measurement times. Defaults to NULL.
#' @param x Vector of measurement values. Defaults to NULL.
#' @param Harg0 If omega_0 is equal to zero (x=0), set omega_0 to Harg0 (w=Harg0). Defaults to NULL.
#' @return A vector with the transformation of the measurement data into the frequency domain.
#' @export
#' @examples This function is called by HMF1compute() and HMFcompute(). This function is not recommended to be used in isolation.
#' Hcompute(omega, time, expr, H0)
#' @references \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005627}
Hcompute <- function(w=NULL,t=NULL,x=NULL,Harg0=NULL) { # continuous Hartley transform approximation
  H <- 0 # initiallized CHT, to be updated in loop over all intervals
  for (i in 1:(length(t)-1)) {
    xf <- x[i+1]
    xi <- x[i]
    tf <- t[i+1]
    ti <- t[i]
    mm <- (xf - xi)/(tf - ti) # slope for linear interpolation
    b <- x[i]-mm*t[i] # intercept for linear interpolation
    
    if (w==0) {w = Harg0}
    Imtcoswtdt <- (mm/(w^2)) * ( (cos(w*tf)+w*tf*sin(w*tf)) - (cos(w*ti)+w*ti*sin(w*ti)) )
    Ibcoswtdt <- (b/w) * ( sin(w*tf) - sin(w*ti) )
    Imtsinwtdt <- (mm/(w^2)) * ( (sin(w*tf)-w*tf*cos(w*tf)) - (sin(w*ti)-w*ti*cos(w*ti)) )
    Ibsinwtdt <- (-b/w) * ( cos(w*tf) - cos(w*ti) )
    H <- H + Imtcoswtdt + Ibcoswtdt + Imtsinwtdt + Ibsinwtdt
  }
  return(H)
}


#' HMF spectral component
#'
#' This function computes the mth HMF spectral component of the measurement profile (Equation 5, Anderson et al., 2017).
#' @param meds Matrix with the first column containing measurement times and the second column containing the corresponding measurements. Defaults to NULL.
#' @param m The m value. 
#' @param W0 The omega_0 term. 
#' @param n Highest order of the system (n=1 for the first degree system that this analysis was designed for).
#' @param Harg0 If omega_0 is equal to zero (x=0), set omega_0 to Harg0 (w=Harg0). 
#' @return The spectral component data.
#' @export
#' @examples This function is not recommended to be used in isolation.
#' HMFcompute(data, m, omega, n, H0)
#' @references \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005627}
HMFcompute <- function(meds,m,W0,n,Harg0) {  # spectral component, i=1
  t <- meds[,1] # time
  x <- meds[,2] # expression data, -ddCt (imputed)
  i <- 1        # order of the derivative, applied for dE/dt (i=1)
  Hhat <- 0     # initialize Hhat, to be updated in loop
  for (j in c(0:n)) {
    casarg <- (i*pi/2) * ((n+m-j)^i)
    cas <- cos(casarg) - sin(casarg)
    Harg <- (n+m-j) * W0 
    H <- Hcompute(Harg,t,x,Harg0)
    Hhat <- Hhat + ((-1)^j) * choose(n,j) * H
  }
  return(Hhat)
} 


#' HMF spectral component of the 1st derivative
#'
#' This function computes HMF spectral component of the 1st derivative of x (dx/dt = x', Equation 6, Anderson et al., 2017).
#' @param meds Matrix with the first column containing measurement times and the second column containing the corresponding measurements. Defaults to NULL.
#' @param m The m value. 
#' @param W0 The omega_0 term. 
#' @param n Highest order of the system (n=1 for the first degree system that this analysis was designed for).
#' @param Harg0 If omega_0 is equal to zero (x=0), set omega_0 to Harg0 (w=Harg0). 
#' @return The spectral component data.
#' @export
#' @examples This function is not recommended to be used in isolation.
#' HMF1compute(data, m, omega, n, H0)
#' @references \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005627}
HMF1compute <- function(meds,m,W0,n,Harg0) {  # spectral component, i=1
  t <- meds[,1] # time
  x <- meds[,2] # expression data, -ddCt (imputed)
  i <- 1        # order of the derivative, applied for dE/dt (i=1)
  Hhat <- 0     # initialize Hhat, to be updated in loop
  for (j in c(0:n)) {
    casarg <- (i*pi/2) 
    cas <- cos(casarg) - sin(casarg)
    Harg <- ((-1)^i) * (n+m-j) * W0 
    H <- Hcompute(Harg,t,x,Harg0)
    Hhat <- Hhat + ((-1)^j) * choose(n,j) * cas * ((n+m-j)^i) * (W0^i) * H
  }
  return(Hhat)
} 


#' X and Y data
#'
#' This function function for getting the X and Y data (Equation 10, Anderson et al., 2017).
#' @param datHMF The data matrix with scaled measurements measurements. Scaled center measurements are provided for 
#' each measurement/annotation (e.g., gene/organ) combination. This format is produced by scale_zeroOne(). Defaults to datHMF.
#' @param Mall A set of m-values. Defaults to NULL.
#' @param Msc A scale factor that determines the M value (M = max(m)) corresponding to the largest range of m-values. 
#' Mmax <- ceiling(Msc * Np / 2) where Np = number of parameters.
#' @param n Highest order of the system (n=1 for the first degree system that this analysis was designed for).
#' @param Harg0 If omega_0 is equal to zero (x=0), set omega_0 to Harg0 (w=Harg0). 
#' @return The X and Y transformations obtained using Equations 4-6 for all m-value ranges
#' @export
#' @examples XY = getXY(datHMF,Mall=Mall,Msc=Msc,n=n,Harg0=Harg0)
#' @references \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005627}
getXY = function(datHMF=datHMF,Mall=NULL,Msc=2,n=1,Harg0=NULL){
  
  W0 = 2*pi / (max(unique(datHMF[,1]))-min(unique(datHMF[,1])))
  
  if(is.null(Mall) == F){Mall = Mall}
  if(is.null(Mall) == T){
    mvals = get_M_vals(datHMF=datHMF,Mall=NULL,Msc=2)
    Mall = mvals$Mall
  }
  
  X <- array( c(0), c(ncol(datHMF)-2, length(unique(datHMF[,2])), length(Mall) ) ) 
  Y <- array( c(0), c(ncol(datHMF)-2, length(unique(datHMF[,2])), length(Mall) ) )
  
  for(anns in 1:length(unique(datHMF[,2]))){ # loop through annotations
    ind_ann = which(datHMF[,2] == unique(datHMF[,2])[anns])
    for(feats in 1:(ncol(datHMF)-2)){ # loop through features
      ind_feat = which(names(datHMF) == names(datHMF)[3:ncol(datHMF)][feats])
      mids = datHMF[ind_ann,c(1,ind_feat)]
      
      # loop for HMF calculations
      for (mr in 1:length(Mall)) {         # loop through every m value
        mm <- Mall[mr]  # this is the m value
        X[feats,anns,mr] <- HMFcompute(mids,mm,W0,n,Harg0)
        Y[feats,anns,mr] <- HMF1compute(mids,mm,W0,n,Harg0)
      } # mr
      
    } # feats, loop through features
  } # anns, loop through annotations
  
  out = list(X=X, Y=Y)
  return(out)
  
} # end getXY()
