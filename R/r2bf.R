
#' Convert correlation \emph{r} to a Bayes Factor
#' 
#' The Bayes Factor is one form of evidence ratio. This function uses default values from the \code{\link[BayesFactor]{ttest.tstat}} function to convert an 
#' effect size \emph{r} to a Bayes factor by first finding the T value, given the sample size \emph{n}. 
#' Can also include as attributes the thresholds, 
#' at that sample size, for conventionally "accepting that there is at least moderate evidence against the tested association" (i.e., BF < .333)
#' or "at least moderate evidence for the tested association" (i.e., BF > 3). Bayes Factors are returned "raw," that is,
#' as the ratio of the evidence for the alternative hypothesis to the evidence for the null hypothesis.
#' All values returned should be treated as \strong{approximate}, and the hypothesis 
#' acceptance thresholds should be treated as \strong{heuristics}. For Bayes Factors for
#' a particular correlation, using \code{\link[BayesFactor]{correlationBF}} is \strong{highly recommended}.
#' 
#' In informal testing, the Bayes Factors estimated using the \code{BayesFactor} package (with default priors) 
#' appear to collapse the 
#' "zone of ambiguity" with sample sizes above around 50 or 75, leading to a large proportion of BF<.333 or BF>3. 
#' This is especially problematic due to the relatively large correlations (e.g., r>.1) identified as "evidence for
#' the null hypothesis." The 
#' user should keep in mind that Bayes Factors [1] are \emph{inherently model comparisons} and [2] are often 
#' \emph{quite dependent on the prior being used} (see, e.g., \code{doi.org/10.1002/wics.1595}), making default,
#' automatic, and interpretable estimation of these evidence ratios a near-impossible task. 
#' 
#' The best alternative
#' is to frame your model comparison in a more formal way, by defining two models (e.g., regressions) with maximally-similar priors
#' and then using bridge sampling to estimate marginal likelihoods with which to calculate evidence ratios (i.e., BFs). 
#' The \code{brms} package includes excellent ways to implement this type of analysis.
#' 
#' For discussion of methods and conventional BF categories, see \code{doi: 10.3758/s13423-012-0295-x}.
#' 
#' @param r 
#' @param n 
#' @param rscale Passes to the \code{rscale} argument of \code{\link[BayesFactor]{ttest.tstat}}. Choose from "medium.narrow", "medium", "wide", and "ultrawide". As in \code{ttest.tstat}, defaults to "medium".
#' @param nsims Not currently used.
#' @param returnThresholds Logical. Should the boundaries indicating "decision thresholds" be estimated?
#'
#' @export
#' 
#' @examples
#' exampleR <- cor(dat_cochraneEtAl_2019_PLOSOne$enum , dat_cochraneEtAl_2019_PLOSOne$mot)
#' r2bf(exampleR , nrow(dat_cochraneEtAl_2019_PLOSOne))
#' 
#' 
#' 
r2bf <- function(r,n,rscale = 'medium',nsims = 20, returnThresholds = FALSE){
  
  require(BayesFactor)
  require(mvtnorm)
  
  getBF <- function(rIn, nSimsIn = nsims){
    
    # this gives much more stable results than the below
    BayesFactor::ttest.tstat(psych::r2t(rIn,n), n, simple = T, rscale = rscale)
    
    ## Getting this to work well is not really happening (estimates for the thresholds are quite high r for both null and non-null inferences, when using default priors)
    # full_bf <- c()
    # for(. in 1:nSimsIn){
    #   
    #   correlDat <- list()
    #   for(curDat in 1:20){
    #     correlDat[[curDat]] <- 
    #       rmvnorm(n,sigma = matrix(c(1,rIn,rIn,1),nrow = 2) )
    #   }
    #   
    #   correlDat <- correlDat[[which.min(abs(unlist(
    #     lapply(correlDat,function(x){cor(x)[1,2]})
    #   )))]]
    #   
    #   full_bf <- c(full_bf
    #                ,correlationBF(correlDat[,1],correlDat[,2])@bayesFactor$bf)
    #   
    # }
    # exp(mean(full_bf))
  }
  
  bf <- getBF(rIn=r)
  
  attr(bf,'results') <- c(n=as.integer(n),r=r,bf = bf)
  
  if(returnThresholds){
    
    nullThresh <- c()
    alternThresh <- c()
    
    for(. in 1:3){
      
      nullThresh <- c(nullThresh
                      ,suppressMessages({
                        suppressWarnings({
                          optimize(function(r_val) {
                            (0.333333333333 - getBF(rIn = r_val,nSimsIn = nsims))^2
                          }, c(1e-06, 0.99999)
                          # ,nsims = nsims
                          )$minimum
                        })
                      })
      )
      
      alternThresh <- c(alternThresh
                        ,suppressMessages({
                          suppressWarnings({
                            optimize(function(r_val,nsims) {
                              (3- getBF(rIn = r_val
                                        # ,nSimsIn = nsims
                                        ))^2
                            }, c(1e-06, 0.99999)
                            # ,nsims = nsims
                            )$minimum
                          })
                        })
      )
      
    }
    
    attr(bf,'results')['nullThresh'] <- median(nullThresh)
    attr(bf,'results')['alternThresh'] <- median(alternThresh)
    
  }
  
  bf
  
}

