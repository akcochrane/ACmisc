#' Apply optimal Yeo-Johnson transformation to a numeric vector
#'
#' Finds the Yeo-Johnson transformation that minimizes the skewness of the vector.
#' Returns the transformed vector with an attribute "lambda" (transformation parameter).
#'
#' @param x numeric vector to transform
#' @param keepMedian Should the result be shifted so that it has the same median as the input vector?
#' @param keepMAD Should the result be scaled so that it has the same dispersion (median absolute deviation) as the input vector?
#'
#' @export
#'
YeoJohn <- function(x,keepMedian=T,keepMAD=T){

  yeojohnson <- VGAM::yeo.johnson

  require(psych)

  xOut <- rep(NA,length(x))
  x_nonNA <- !is.na(x)

  x <- na.omit(x)

  x_expanded <- quantile(x,seq(.0001,.9999,length = 1E4))

  med_orig <- median(x_expanded)
  mad_orig <- mad(x_expanded)
  if(mad_orig == 0){mad_orig <- sd(x_expanded)}

  rootFun <- function(lambda){return(skew(yeojohnson(x_expanded,lambda))^2)}
  bestLambda <- optimize(rootFun,interval = c(-20,20))

  xOut[x_nonNA] <- yeojohnson(x,bestLambda$minimum)
  lambda <- bestLambda$minimum
  rawDat <- x
  rawDat_skew <- skew(x)
  transformedDat <- xOut

  if(keepMedian || keepMAD){

    med_transf <- median(xOut,na.rm=T)

    xOut <- xOut - med_transf

    if(keepMAD){
      mad_new <- mad(xOut,na.rm=T)
      if(mad_new == 0){mad_new <- sd(xOut,na.rm=T)}
      xOut <- xOut/mad_new
      xOut <- xOut*mad_orig}

    if(keepMedian){
      xOut <- xOut+med_orig
    }else{xOut <- xOut + med_transf}


    attr(xOut,'transformedDat') <- transformedDat
  }

  attr(xOut,'lambda') <- lambda
  attr(xOut,'rawDat') <- rawDat
  attr(xOut,'rawDat_skew') <- rawDat_skew

  return(xOut)

  if(F){
    x <- rnorm(200,10,2)
    skewYJ <- data.frame() ; for(curPow in seq(-4,4,.01)){
      curYJ <- YeoJohn(x^curPow)
      skewYJ <- rbind(skewYJ,data.frame(pow = curPow,skew = attr(curYJ,'rawDat_skew'),lambda = attr(curYJ,'lambda')))
      }
    psych::pairs.panels(skewYJ)
    }
}
