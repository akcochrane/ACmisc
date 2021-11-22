#' Normalize a vector to be between .001 and .999
#' 
#' To force a variable to have an easily interpretable min and max.
#' 
#' @param dat numeric vector to normalize
#' @keywords  keywords
#' @export
#' @examples 
#' examples
#' zeroOneNorm(rnorm(30))
#' zeroOneNorm(dat_cochraneEtAl_2019_PLOSOne$mot)

zeroOneNorm <- function(dat){
  
  # dat <- rnorm(30,0,4)
  
  dOut <- dat - min(dat,na.rm=T)
  dOut <- dOut/max(dOut,na.rm=T)
  dOut <- dOut/(1/.998)
  dOut <- dOut + .001
  
  return(dOut)
}


