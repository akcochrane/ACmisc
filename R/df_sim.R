#' Simulate new data from a numeric and/or logical data frame or matrix
#'
#' Simulates a new dataset, conforming to a reference dataset as closely as possible
#' By resampling the existing data, only possible values are allowed in simulating
#' a new dataset. This function attempts to minimize error in correlational
#' structure, median, and MAD.
#'
#' Error, that is minimized, is defined as:
#' \code{err <-
#' sum((orig_correl - cur_correl)^2)*mean(orig_mad) +
#' sum((orig_med - cur_med)^2)+
#' sum((orig_mad - cur_mad)^2)}
#'
#' Where \code{orig_correl} and \code{cur_correl} are the original and current correlation matrices.
#' This is a somewhat arbitrary loss function, but it should do OK for matching new and
#' old datasets.
#'
#' Next step: take a groupingVar argument, and simulate from subsets of the data for each
#' unique value of groupingVar, then concatenate those subsets of data (while generating
#' random identifiers for each group? Probably should have an anonymize=T argument as well)
#'
#' @param nsim Number of cases to simulate
#' @param x Input data frame or matrix. Must be numeric or logical
#' @param cor_method Minimize error in 'spearman' or 'pearson' correlations?
#' @param nTries This number of resamples are iteratively proposed, and the one that best matches the original is retained.
#'
#' @export
#'
#' @examples
#' df_raw <- data.frame(x = rnorm(20),y=rep(c(1,2,3,4),5))
#' df_raw$z <- round((df_raw$x + df_raw$y + rnorm(20,-3,3))*4)/4
#'
#' df_sim <- sim_dat(20,df_raw)
#' abs(cor(df_raw) - cor(df_sim))
#' apply(df_raw,2,mean) ; apply(df_sim,2,mean)
#' apply(df_raw,2,sd) ; apply(df_sim,2,sd)
#'
sim_dat <- function(nsim,x,cor_method='spearman',nTries = 200){

  ## ## get initial sample with decently matched correlational structure
{
  orig_correl <- suppressWarnings({cor(x,method=cor_method,use='pairwise')})
  orig_med <- apply(x,2,median,na.rm=T)
  orig_mad <- apply(x,2,mad,na.rm=T)


  best_err <- 1E15 ;  for(curStart in 1:nTries){ # 1st through nTriesth initialization attempts
    curSim <- apply(x,2,sample,10,replace=T)

    cur_correl <- suppressWarnings({cor(curSim,method=cor_method,use='pairwise')})
    cur_correlErr <- sum((orig_correl - cur_correl)^2)

    cur_med <- apply(curSim,2,median,na.rm=T)
    cur_mad <- apply(curSim,2,mad,na.rm=T)

    cur_err <-
      cur_correlErr*mean(orig_mad) +
      sum((orig_med - cur_med)^2)+
      sum((orig_mad - cur_mad)^2) +
      1E15 * (nrow(na.omit(curSim)) < nrow(na.omit(x[sample(nrow(x),nrow(curSim),replace=T),])))# ensure NAs don't propagate

    if(cur_err < best_err){
      bestSim <- curSim
      best_err <- cur_err
    }
  }
}

  while(nrow(bestSim) < nsim){

    best_err <- 1E15 ; for(. in 1:nTries){
    curSim <- rbind(bestSim,apply(x,2,sample,1,replace=T))
    cur_correl <- suppressWarnings({cor(curSim,method=cor_method,use='pairwise')})
    cur_correlErr <- sum((orig_correl - cur_correl)^2)

    cur_med <- apply(curSim,2,median,na.rm=T)
    cur_mad <- apply(curSim,2,mad,na.rm=T)

    cur_err <-
      cur_correlErr*mean(orig_mad) +
      sum((orig_med - cur_med)^2)+
      sum((orig_mad - cur_mad)^2) +
      1E15 * (nrow(na.omit(curSim)) < nrow(na.omit(x[sample(nrow(x),nrow(curSim),replace=T),]))) # ensure NAs don't propagate


    if(cur_err < best_err){
      bestTmpSim <- curSim
      best_err <- cur_err
    }
    }
    bestSim <- bestTmpSim ; rm(bestTmpSim)
  }

  ## look at performance:
  # apply(bestSim,2,mean,na.rm=T) ; apply(x,2,mean,na.rm=T)
  # apply(bestSim,2,sd,na.rm=T) ; apply(x,2,sd,na.rm=T)
  # cor(bestSim) ; cor(x)

  return(bestSim[1:nsim,])
}
