
#' Plot simulation data
#'
#' This function plots simulation data
#' @param negddCt Raw dara file.
#' @param Sims Simulation result.
#' @export
#' @examples 
#' sim_plots(multiorgan, result$best_sim)
sim_plots = function(negddCt=NULL,Sims=NULL){
  
  pltInd = 2
  Nr = length(unique(negddCt[,2]))
  Ng = length(ncol(negddCt)[3:ncol(negddCt)])
  for (r in 1:Nr) { # organ loop, r
    for (g in 1:Ng){ # gene loog, g
      
      # isolate experimental data
      indr = which(negddCt[,2]==unique(negddCt[,2])[r])
      indg = which(colnames(negddCt) == colnames(negddCt)[3:ncol(negddCt)][g])
      dat = negddCt[indr,c(1,2,indg)]
      
      # compute data for plotting (mean and standard error)
      err <- matrix(c(0),length(unique(dat[,1])),3)
      for(ii in 1:nrow(err)){
        ind_age = which(dat[,1] == unique(dat[,1])[ii])
        err[ii,1] = unique(dat[,1])[ii] # age
        err[ii,2] = mean(dat[ind_age,3]) # mean
        err[ii,3] = sd(dat[ind_age,3]) / sqrt(length(dat[ind_age,3]))
      }
      
      ### plot data and simulation
      x <- err[,1]; y <- err[,2]; sem <- err[,3]
      sim <- Sims[,pltInd]
      OldMax=1; OldMin=0; NewMax=max(y); NewMin=min(y)
      OldRange = (OldMax - OldMin); NewRange = (NewMax - NewMin)  
      simScaled = (((sim - OldMin) * NewRange) / OldRange) + NewMin
      sc = (max(dat[,3]) - min(dat[,3]))/20
      plot(dat[,1],dat[,3],xlab="Age (wk)",col="red",ylab=paste(names(dat)[3],"expression (-ddCt)" ),
           ylim=c((min(dat[,3])-sc),(max(dat[,3])+sc)),
           xlim=c(min(dat[,1])-1,max(dat[,1])+1),cex.lab=1.2,cex.main=1.5 )
      title( paste(dat[1,2]), line = -1)
      for (ii in 1:length(x)) {
        if (sem[ii]>0) {
          Hmisc::errbar(x[ii], y[ii], y[ii]+sem[ii], y[ii]-sem[ii], add=T, pch=1, cex=2, cap=.05, col="black", errbar.col="black")
        }
      }
      lines(Sims[,1],simScaled,col="black")
      
      pltInd = pltInd + 1
      
    } # gene loog, g
  } # organ loop, r
  
}