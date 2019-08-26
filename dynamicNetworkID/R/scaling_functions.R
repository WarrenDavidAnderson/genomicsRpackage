
#' Normalization 0,1
#'
#' This function normalizes input data to the range (0,1).
#' @param mat Matrix with the first column containing measurement times and the second column 
#' containing annotations (organs). The third and subsequent columns contain measurement data.
#' @param center The metric used to compute the center, options include "mean" and "median" 
#' @return The data matrix with scaled measurements measurements. Scaled center measurements are provided for 
#' each measurement/annotation (e.g., gene/organ) combination. This is the data format required for HMF().
#' @export
#' @examples datHMF = scale_zeroOne(matrix, center="mean")
scale_zeroOne <- function(mat,center="mean") {
  
  # loop through annotations
  all_centers = c()
  annots = unique(mat[,2])
  for(ii in 1:length(annots)) {
    ind_ann = which(mat[,2] == annots[ii])
    x0 = mat[ind_ann,] # data set for a specific annotation
    time = x0[,1]   # this function assumes that column 1 of x0 is time
    xscaled = x0    # this function assumes that column 3 and on are dynamic variables
    centers = matrix(0, length(unique(time)), ncol(x0)) 
    centers = as.data.frame(centers)
    centers[,1] = unique(time)
    centers[,2] = unique(x0[,2])
    names(centers) = names(x0)
    for(i in 1:length(unique(time))){ # loop through each time and compute the center
      indT = which(x0[,1] == unique(time)[i]) # indices for time i
      if(center=="mean") meanT = apply(x0[indT,3:ncol(x0)], 2, mean)      # vector of means, one for each feature, at time i
      if(center=="median") meanT = apply(x0[indT,3:ncol(x0)], 2, median)  # vector of medians, one for each feature, at time i
      centers[i,3:ncol(x0)] = meanT
    } # i, loop through times
    minT = apply(centers[,3:ncol(x0)], 2, min) # vector of minimums (each feature) at time i
    maxT = apply(centers[,3:ncol(x0)], 2, max) # vector of maximums (each feature) at time i
    mins = do.call("rbind", replicate(nrow(centers), minT, simplify=F))  # matrix of minimums
    maxs = do.call("rbind", replicate(nrow(centers), maxT, simplify=F))  # matrix of maximums
    xscaled = centers
    xscaled[,3:ncol(centers)] = (centers[,3:ncol(centers)] - mins) / (maxs - mins) # apply 0,1 scaling
    all_centers = rbind(all_centers, xscaled)
  } # ii, annotation loop
  return(all_centers)
  
}