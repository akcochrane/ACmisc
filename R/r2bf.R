
#' Convert correlation \emph{r} to a Bayes Factor
#' 
#' Uses default values from the \code{BayesFactor::ttest.tstat} function to convert an 
#' effect size \emph{r} to a Bayes factor, given the sample size \emph{n}. Also includes as attributes the thresholds, 
#' at that sample size, for conventionally "accepting the null hypothesis" (i.e., BF < .333)
#' or "rejecting the null hypothesis" (i.e., BF > 3). Bayes Factors are returned "raw," that is,
#' as the ratio of the evidence for the alternative hypothesis to the evidence for the null hypothesis.
#'
#' @param r 
#' @param n 
#'
#' @export
#'
#' @examples
#' exampleR <- cor(dat_cochraneEtAl_2019_PLOSOne$enum , dat_cochraneEtAl_2019_PLOSOne$mot)
#' r2bf(exampleR , nrow(dat_cochraneEtAl_2019_PLOSOne))
#' 
r2bf <- function(r,n){
  
  getBF <- function(r){
    BayesFactor::ttest.tstat(psych::r2t(r, 
                                        n), n, simple = T)
  }
  
 bf <-  getBF(r=r)
 
 attr(bf,'results') <- c(n=as.integer(n),r=r,bf = bf)
 
 attr(bf,'results')['nullThresh'] <- suppressMessages({
    suppressWarnings({
      optimise(function(r_val) {
        (0.333333333333 - getBF(r_val))^2
      }, c(1e-06, 0.99999))$minimum
    })
  })
 
 attr(bf,'results')['alternThresh'] <- suppressMessages({
   suppressWarnings({
     optimise(function(r_val) {
       (3- getBF(r_val))^2
     }, c(1e-06, 0.99999))$minimum
   })
 })

   
  
 
 bf
  
}

