
#' Find by-trial Wiener Diffusion Model parameters
#' 
#' Estimates a set of Wiener Diffusion Model parameters for a vector of formatted
#' data. Estimation uses \code{\link{RWiener::wdm}} and data must be formatted as a 
#' 1-dimensional vector according the to \code{RWiener} package (i.e., RT * 2 
#' * (boundary=='upper' - 0.5) ). 
#'
#' @param dat Response time vector. Positive values indicates upper-boundary while negative values indicate lower-boundary.
#' @param fit_DR Logical or numeric scalar. If \code{TRUE} then drift rate is estimated for each trial, if \code{FALSE} only the overall parameter is estimated, and if numeric then the parameter is fixed to that value.
#' @param fit_BS Logical or numeric scalar. If \code{TRUE} then boundary separation is estimated for each trial, if \code{FALSE} only the overall parameter is estimated, and if numeric then the parameter is fixed to that value.
#' @param fit_NDT Logical or numeric scalar. If \code{TRUE} then non-decision time is estimated for each trial, if \code{FALSE} only the overall parameter is estimated, and if numeric then the parameter is fixed to that value.
#' @param fit_Bias Logical or numeric scalar. If \code{TRUE} then bias is estimated for each trial, if \code{FALSE} only the overall parameter is estimated, and if numeric then the parameter is fixed to that value.
#' @param winsorize Optional numeric scalar. If supplied, then all univariate values greater than number of MAD from the median are changed to be instead that number of MAD from the median (using \code{mad()} defaults).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- (rnorm(200,.2,.05) + rexp(200,1 / .2))*sample(c(-1,1,1),200,replace=T)
#'     m <- trialWDM(dat)
#'     
#'     m2 <- trialWDM(dat, fit_NDT = F , winsorize = 2)
#'     
#'     }
trialWDM <- function(dat,fit_DR = T , fit_BS=T , fit_NDT = T , fit_Bias = F , winsorize = NA){


require(RWiener)

overallMod <-  wdm(dat 
                   , alpha = if(is.logical(fit_BS)){NULL}else{fit_BS ; fit_BS <- F}
                   ,beta = if(is.logical(fit_Bias)){NULL}else{fit_Bias ; fit_Bias <- F}
                   ,tau = if(is.logical(fit_NDT)){NULL}else{fit_NDT ; fit_NDT <- F}
                   ,delta = if(is.logical(fit_DR)){NULL}else{fit_DR ; fit_DR <- F}
)

overallPars <- overallMod$coefficients

drop1length <- length(dat)-1

# # for reference:
# trialChange = (big-small) * length(small)
# trialspecific = (big-small)*length(small) + big
# # demo, for the mean:
# x <- 1:10 ; x_m <- mean(x) ; x_recovered <- c() ;  for(curDrop in x){curVect <- x[-curDrop] ; curChange <- x_m - mean(curVect) ; x_recovered <- c(x_recovered, x_m + curChange * length(curVect) )}


trialWDM <- data.frame() ; sddf <- data.frame() ; for(curRow in 1:length(dat)){
  
curWDM <- wdm(dat 
              , alpha = if(fit_BS){NULL}else{overallPars['alpha']}
              ,beta = if(fit_Bias){NULL}else{overallPars['beta']}
              ,tau = if(fit_NDT){NULL}else{overallPars['tau']}
              ,delta = if(fit_DR){NULL}else{overallPars['delta']}
)
  
  trialWDM <- rbind(trialWDM,data.frame(
    DR = drop1length*(overallPars['delta'] - curWDM$coefficients['delta']) + overallPars['delta']
    ,BS =  drop1length*(overallPars['alpha'] - curWDM$coefficients['alpha']) + overallPars['alpha']
    ,NDT = drop1length*(overallPars['tau'] - curWDM$coefficients['tau'] ) + overallPars['tau'] 
    ,bias = drop1length*(overallPars['beta'] - curWDM$coefficients['beta'] ) +overallPars['beta']
  ))
  
  # sddf <- rbind(sddf, t(data.frame(summary(curWDM)$sd)))
  
  rm(curWDM)
  
}

# SDs <- apply(sddf, 2, mean, trim = .25) # hopefully more robust than the overall point estimate
# 
# for(curPar in names(SDs)){
#   
#   curParName <- switch (curPar,
#     'alpha' = 'BS'
#     ,'beta' = 'bias'
#     ,'delta' = 'DR'
#     ,'tau' = 'NDT'
#   )
#  
#   trialWDM[,curParName] <-  # rescale to the "true" SD 
#     scale(trialWDM[,curParName])*SDs[curPar] + mean(trialWDM[,curParName])
#   
# }

if(is.numeric(winsorize)){

  trialWDM <- apply(trialWDM,2,FUN=function(x){
  xmed <- median(x)
  xmad <- mad(x)
  
  x[x > (xmed + winsorize*xmad)] <- (xmed + winsorize*xmad)
  x[x < (xmed - winsorize*xmad)] <- (xmed - winsorize*xmad)
  # x
})
}

trialWDM <- data.frame(raw_data = dat[1:nrow(trialWDM)], trialWDM)

attr(trialWDM, 'overall_model') <- overallMod

return(trialWDM)
}
