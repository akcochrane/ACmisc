
#' Convert correlation \emph{r} to a Bayes Factor
#' 
#' Uses default values from the \code{BayesFactor::correlationBF} function to convert an 
#' effect size \emph{r} to a Bayes factor by simulating datasets, given the sample size \emph{n}. Can also include as attributes the thresholds, 
#' at that sample size, for conventionally "accepting the null hypothesis" (i.e., BF < .333)
#' or "rejecting the null hypothesis" (i.e., BF > 3). Bayes Factors are returned "raw," that is,
#' as the ratio of the evidence for the alternative hypothesis to the evidence for the null hypothesis.
#' All values returned should be treated as \strong{approximate}, and the hypothesis 
#' acceptance thresholds should be treated as \strong{heuristics}. For Bayes Factors for
#' a particular correlation, using \code{\link[BayesFactor]{correlationBF}} is \strong{highly recommended}.
#' 
#' @param r 
#' @param n 
#' @param rscale Passes to the \code{rscale} argument of \code{\link[BayesFactor]{correlationBF}}. Choose from "medium.narrow", "medium", "wide", and "ultrawide". As in \code{correlationBF}, defaults to "medium".
#' @param nsims Number of simulated datasets for which to find the Bayes Factors; more simulations will lead to more-stable results
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
    
    ## previously this method, but it seems to take a shortcut (but maybe gives the same result in the end..?)
    # BayesFactor::ttest.tstat(psych::r2t(r,n), n, simple = T)
    
    full_bf <- c()
    for(. in 1:nSimsIn){
      
      correlDat <- list()
      for(curDat in 1:20){
        correlDat[[curDat]] <- 
          rmvnorm(n,sigma = matrix(c(1,rIn,rIn,1),nrow = 2) )
      }
      
      correlDat <- correlDat[[which.min(abs(unlist(
        lapply(correlDat,function(x){cor(x)[1,2]})
      )))]]
      
      full_bf <- c(full_bf
                   ,correlationBF(correlDat[,1],correlDat[,2])@bayesFactor$bf)
      
    }
    exp(mean(full_bf))
  }
  
  bf <- getBF(rIn=r)
  
  attr(bf,'results') <- c(n=as.integer(n),r=r,bf = bf)
  
  if(returnThresholds){
    
    nullThresh <- c()
    alternThresh <- c()
    
    for(. in 1:25){
      
      nullThresh <- c(nullThresh
                      ,suppressMessages({
                        suppressWarnings({
                          optimize(function(r_val,nsims) {
                            (0.333333333333 - getBF(rIn = r_val,nSimsIn = nsims))^2
                          }, c(1e-06, 0.99999)
                          ,nsims = nsims)$minimum
                        })
                      })
      )
      
      alternThresh <- c(alternThresh
                        ,suppressMessages({
                          suppressWarnings({
                            optimize(function(r_val,nsims) {
                              (3- getBF(rIn = r_val,nSimsIn = nsims))^2
                            }, c(1e-06, 0.99999)
                            ,nsims = nsims)$minimum
                          })
                        })
      )
      
    }
    
    attr(bf,'results')['nullThresh'] <- median(nullThresh)
    attr(bf,'results')['alternThresh'] <- median(alternThresh)
    
  }
  
  bf
  
}

