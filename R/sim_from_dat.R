
#' Simulate new data from a dataset's quantiles
#'
#' @param nSim Number of simulated observations
#' @param datIn Data to simulate from
#' @param minPerClust [Optional] If subsets of the dataset are desired to be simulated from, these subsets can be identified through K-means clustering. See Details.
#' @param nTries Number of simulated datasets to generate. Only one (i.e., the one with the rank correlation closest to the original dataset's) will be returned.
#'
#' @details 
#' 
#' Given an existing dataset's rank correlation structure and the quantiles of 
#' the variables in that dataset, a new dataset is simulated from the quantiles
#' and attempting to match the rank correlation of the original.
#' 
#' Sometimes it may be inappropriate to assume that rank correlations would be 
#' homogeneous across an entire dataset, and instead there may be subsets of the 
#' data that show different patterns than the full dataset (e.g., Simpson's paradox).
#' If \code{minPerClust} is less than the number of observations in \code{datIn},
#' then k-means clustering (default R \code{\link{kmeans}}) is used iteratively to 
#' determine the maximum number of clusters for which the minimum cluster size
#' is at least \code{minPerClust}. If there are at least two clusters satisfying
#' this criterion, then the simulation from quantiles and rank correlations is 
#' completed for each identified cluster separately, and these are concatenated in the 
#' returned data frame. In general, the smaller the value of \code{minPerClust},
#' [1] the closer the simulated data will be to the original data (including 
#' undesirable noise or other idiosyncrasies), and [2] the longer the run time
#' will be. Given that correlations are unlikely to be reliable with small numbers
#' of observations, it would be very strange to have \code{minPerClust} be below
#' 20.
#' 
#' Requires the \code{mvtnorm} package.
#'
#' @export
#' 
#' @seealso 
#' \code{\link{sim_dat}} accomplishes something similar, although that function resamples directly from the
#' original data while this one resamples from the quantiles of the original data.
#'
#' @examples
#' 
#' # note that non-numeric variables will be dropped
#' new_iris_1 <- sim_from_dat(250,iris) 
#' 
#' # For the iris dataset, this will simulate from separate empirically-defined
#' # clusters to better reflect the clustered nature of the original data.
#' new_iris_2 <- sim_from_dat(250,iris, minPerClust = 40) 
#' 
sim_from_dat <- function(nSim, datIn, minPerClust = 100, nTries = 100){
  
  library(mvtnorm)
  
  # datIn <- iris
  
  datIn <- suppressWarnings(datIn[,!is.na(apply(datIn,2,sd))])
  
  ## find the largest number of clusters, for which observations per cluster are more than minPerClust
  if(nrow(datIn) > minPerClust){
    minClust <- 1E20 ; nClust <- 1 ; clustList <- list()
    while(minClust > minPerClust){
      
      clustList[[nClust]] <- kmeans(datIn, nClust, nstart = 10)
      minClust <- min(clustList[[nClust]]$size)
      nClust <- nClust + 1
    }
    datClusts <- clustList[[length(clustList) - 1]]$cluster
  }else{
    datClusts <- rep(1,nrow(datIn))
  }
  
  nPerClust <- round((table(datClusts)/sum(table(datClusts)) ) * nSim)
  while(sum(nPerClust) < nSim){ nPerClust[1] <- nPerClust[1] + 1 } # deal with rounding error
  
  ## get a simulated dataset for each cluster, and concatenate
  datOut <- data.frame()
  for(curClust in unique(datClusts)){
    
    datIn_clust <- datIn[datClusts == curClust,]
    
    ## find the original Spearman correlations
    orig_spearman <- cor(datIn_clust,use= 'pairwise',method='spearman')
    
    tries <- replicate(nTries,{
      
      new_mv_quant <- data.frame(pnorm(rmvnorm(nPerClust[curClust], sigma = orig_spearman)))
      colnames(new_mv_quant) <- colnames(datIn_clust)
      
      new_mv_list <- list() ; for(curVar in colnames(datIn_clust)){
        new_mv_list[[curVar]] <- quantile(datIn_clust[,curVar],new_mv_quant[,curVar])
      }
      
      do.call('cbind',new_mv_list)
    },simplify = F)
    
    corDeviance <- unlist(
      lapply(tries,function(x){sum((orig_spearman - cor(x,method= 'spearman'))^2)})
    )
    datOut <- rbind(datOut
                    ,data.frame(tries[[which.min(corDeviance)]]) 
    )
    rm(datIn_clust, orig_spearman, tries,corDeviance)
  }
  
  rownames(datOut) <- NULL
  
  attr(datOut, 'nClusters') <- length(unique(datClusts))
  attr(datOut, 'cluster_id_original_dat') <- datClusts
  
  return(datOut) 
}
