
#' Simulate new data from a dataset's quantiles
#'
#' @param nSim Number of simulated observations
#' @param datIn Data to simulate from
#' @param nTries Number of simulated datasets to generate. Only one, the one with the rank correlation closest to the original dataset's, will be returned.
#'
#' Given an existing dataset's rank correlation structure and the quantiles of 
#' the variables in that dataset, a new dataset is simulated from the quantiles
#' and attempting to match the rank correlation of the original.
#'
#' @export
#' 
#' @seealso 
#' \code{\link{sim_dat}} accomplishes something similar, although that function resamples directly from the
#' original data while this one resamples from the quantiles of the original data.
#'
#' @examples
#' 
#' new_iris <- sim_from_dat(200,iris) # note that non-numeric variables will be dropped
#' 
sim_from_dat <- function(nSim, datIn, nTries = 500){
  
  library(mvtnorm)
  
  # datIn <- iris
  
  datIn <- suppressWarnings(datIn[,!is.na(apply(datIn,2,sd))])
  
  orig_spearman <- cor(datIn,use= 'pairwise',method='spearman')
  
  tries <- replicate(nTries,{
    
    new_mv_quant <- data.frame(pnorm(rmvnorm(nSim, sigma = orig_spearman)))
    colnames(new_mv_quant) <- colnames(datIn)
    
    new_mv_list <- list() ; for(curVar in colnames(datIn)){
      new_mv_list[[curVar]] <- quantile(datIn[,curVar],new_mv_quant[,curVar])
    }
    
    do.call('cbind',new_mv_list)
  },simplify = F)
  
  corDeviance <- unlist(
    lapply(tries,function(x){sum((orig_spearman - cor(x,method= 'spearman'))^2)})
  )
  
  datOut <-  data.frame(tries[[which.min(corDeviance)]])
  rownames(datOut) <- NULL
  
  return(datOut) 
}
