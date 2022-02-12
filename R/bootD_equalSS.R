#' Bootstrap equal-sample-size Cohen's d
#' 
#' Unequal sample sizes can skew estimates of group differences. This function
#' takes two numeric vectors, \code{x} and \code{y}, and repeatedly calculates 
#' Cohen's d for the difference in means. On each iteration the smaller of the two vectors 
#' is used whole, while the larger of the two vectors is randomly subsampled with the same number
#' of elements as the smaller vector (NAs excluded). Several 
#' quantiles of the resampled distributions of Cohen's d
#' are returned, along with the true mean difference.
#'
#' @param x The first vector
#' @param y The second vector
#' @param nBoot The number of bootstrap iterations
#'
#' @export
#'
#' @examples
#' x <- c(NA,rnorm(3000,.5))
#' y <- rnorm(30)
#' bootD_equalSS(x,y)
#' 
bootD_equalSS <- function(x,y,nBoot = 2000){
  library(psych)
  x <- na.omit(x)
  y <- na.omit(y)
  
  lengthX <- length(x)
  lengthY <- length(y)
  
  if(lengthX > lengthY){
    longerVar <- x
    shorterVar <- y
  }else{
    longerVar <- y
    shorterVar <- x
  }
  
  sampleN <- length(shorterVar)
  
  cohenDs <- replicate(5,
                       replicate(ceiling(nBoot/5),{
                         curDat <- rbind(
                           data.frame(curSample = shorterVar, groupVar  = 0)
                           ,data.frame(curSample = sample(longerVar,sampleN), groupVar  = 1)
                         )
                         
                         ## Could calculate d by hand faster, but I want to make sure the 
                         ## exact correct formula (same as psych package) is being used
                         
                         # cohen.d(curDat$curSample,curDat$groupVar)$cohen.d[1,'effect']
                         
                         ## this is faster than the above, and provides the same results
                         r2d(cor(curDat$curSample,curDat$groupVar))
                         
                       },simplify = 'array')
  )
  
  trueMeanDiff <- mean(x,na.rm=T) - mean(y,na.rm=T)
  
  if(sign(trueMeanDiff) != sign(mean(cohenDs)) ){
    cohenDs <- -cohenDs}
  
  cohenD_quants <- apply(cohenDs,2,FUN = quantile, probs = c(.005,.025,.25,.5,.75,.975,.995))
  

  
  cohenD_output <- data.frame(trueMeanDiff = trueMeanDiff,
                              nPerGroup = sampleN,
                              t(apply(cohenD_quants,1,mean))
  )
  colnames(cohenD_output) <- gsub('X','Q',colnames(cohenD_output))
  
  return(cohenD_output)
}


