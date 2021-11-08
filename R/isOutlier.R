

#' Find multivariate outliers (Leys et al., 2018, JESP)
#'
#' Tests for multivariate outliers using a robust method (cited above). In summary, a minimum proportion of
#' cases are used to estimate the location and scatter of the multivariate space, then that location and scatter are
#' used to test for [multivariate] outlying values.
#'
#' Takes a numeric data frame. Returns a logical vector indicating whether the cases should be excluded.
#'
#' Is \emph{NOT} robust to missing data.
#'
#' @seealso
#'  \code{\link[MASS]{cov.rob}} for estimating the robust location and scatter.
#'
#' @param dat data
#' @param dat_proportion Proportion of data to retain in estimating the robust multivariate location and scatter. Defaults to testing at various proportions.
#' @param reject_alpha P-value threshold for excluding cases (using a chi-square test on Mahalanobis distances; Leys et al, 2018)
#'
#' @export
#'
isOutlier <- function(dat,dat_proportion = 'check',reject_alpha = .01){


  library(MASS)


  findOutliers <- function(dat,dat_proportion,reject_alpha=reject_alpha){
    covar_fit <- cov.mcd(dat,quantile.used = nrow(dat)*dat_proportion,nsamp='best')  # covariances of at least ##% of the data
    mhmcd <-mahalanobis(dat,covar_fit$center,covar_fit$cov)

    cutoff <- (qchisq(p = 1-reject_alpha, df = ncol(dat)))

    excluded <- rep(F,length(mhmcd))

    excluded[ which( mhmcd > cutoff ) ] <- T

    names(excluded) <- rownames(dat)

    excluded
  }


  if(is.numeric(dat_proportion)){return(findOutliers(dat,dat_proportion = dat_proportion,reject_alpha))
  }else{

    excluded <- data.frame(prop_70 = rep(NA,nrow(dat)))

    for(curProp in seq(.70,.98,.02)){

      curOutliers <- rep(NA,nrow(dat))

      try({
        curOutliers <- findOutliers(dat,curProp,reject_alpha)
      },silent=T)

      excluded[,paste0('prop_',curProp*100)] <- curOutliers

      rownames(excluded) <- rownames(dat)

    }

    comment(excluded) <- 'Columns represent proportions of the data retained in estimating robust covariance matrices
  (i.e., the centroid and spread). See MASS::cov.mcd.'

      attr(excluded,'tradeoff') <- data.frame(estim_prop = seq(.7,.98,.02),excl_prop = colMeans(excluded))

    return(excluded)
  }


}

