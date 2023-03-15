
#' Fit multivariate skew-T distribution to a numeric dataset
#'
#' @param d Input data. Must be numeric, and will likely fail if there are fewer than three unique values in any column. In general, the function is designed for data with plausibly skew-T distributions.
#' @param nRanGen Confidence intervals are estimated using random variates of the multivariate skew-T distribution. This argument defines how many random samples to generate in the estimation.
#'
#' This function treats a dataset as being multivariate skew-T distributed (see \code{\link[sn]{dmst}}),
#' and fits this distribution. Given the distributional assumption, various useful 
#' things can be used; these include the central tendency and CI (the main output 
#' of the function; central tendency is the mean of the middle 1 percent of generated
#' random variates, calculated using \code{mean(., trim = .495)}). Output also includes 
#' various \code{attributes}, including the
#' correlation matrix (\code{rob_correl}), or the 
#' actual multivariate skew-T parameters (\code{mvt_coef}).
#' 
#' Because distributional quantiles (CI) are estimated using randomly-generated
#' values from the multivariate distribution, increasing the number of values
#' will make the estimates more stable. Conversely, smaller values of \code{nRanGen}
#' will make estimates less stable.
#'
#' Relies on the \code{\link[sn]{`sn-package`}}) package. 
#'
#' @export
#'
#' @examples
#' dat <- data.frame(x = rexp(200) , y= exp(rnorm(200)) , z = log(rnorm(200 , 5)))
#' got_mvt <- dat2mvt(dat)
#' got_mvt ; apply(dat, 2, quantile , c(.025 , .5 , .975))
#' 
dat2mvt <- function(d, nRanGen = 5E4){
  
  # what would be good would be to have two types of example:
  # 1] showing that it can indeed recover the correlations and skews, when generating from the true distribution
  # 2]a-c] showing that it's better than assuming Gaussian-ness, when outliers, censoring, or skews are present
  
  library(sn)
  
  mvt_fit <- selm(as.formula(paste0('cbind(',paste(colnames(d),collapse=','),') ~ 1'))
                  ,data = d
                  ,family = 'ST')  
  
  # may fail:
  suppressMessages({
    mvt_coef <- coef(mvt_fit , param.type = 'CP', vector = F)
  })
  if(is.null(mvt_coef)){
    mvt_coef <- coef(mvt_fit , param.type = 'pseudo-CP', vector = F) # will work
    names(mvt_coef) <- gsub('~','',names(mvt_coef),fixed = T)
  }
  
  dpM <- list(xi = as.numeric(mvt_fit@param$dp$beta)
              , Omega= mvt_fit@param$dp$Omega
              , alpha= mvt_fit@param$dp$alpha
              , nu=  mvt_fit@param$dp$nu # only for mvt
  )
  
  ranGen <- rmst(nRanGen, dp=dpM)
  
  meanMiddle1percent_CI95 <- data.frame(
    quant025 = apply(ranGen, 2, quantile, .025)
    ,meanMiddle1percent = apply(ranGen, 2, mean , trim = .495)
    ,quant975 = apply(ranGen, 2, quantile , .975)
  )
  rownames(meanMiddle1percent_CI95) <- colnames(d)
  
  attr(meanMiddle1percent_CI95 , 'rob_correl') <- 
    cov2cor(mvt_coef$var.cov)
  
  try({
    attr(meanMiddle1percent_CI95 , 'rob_correl_pVal') <- 
      psych::r.test(
        mvt_fit@size['n.obs'] - (2 + mvt_fit@size['n.param'])
        ,attr(meanMiddle1percent_CI95 , 'rob_correl')
      )$p
  },silent = T)
  
  attr(meanMiddle1percent_CI95 , 'mvt_coef') <-   
    mvt_coef
  
  attr(meanMiddle1percent_CI95 , 'mvt_fit') <-   
    mvt_fit
  
  attr(meanMiddle1percent_CI95 , 'mvt_dp') <-   
    dpM
  
  attr(meanMiddle1percent_CI95 , 'mvt_sim_instructions') <- 
    "use the function sn::rmst(n,dp=attr(output, 'mvt_dp'))"
  
  return(meanMiddle1percent_CI95 )
  
}
