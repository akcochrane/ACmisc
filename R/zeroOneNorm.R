

#' Normalize a vector to be between .001 and .999
#' 
#' Description 
#' 
#' Details
#' 
#' @param dat numeric vector to normalize
#' @keywords  keywords
#' @export
#' @examples 
#' example

zeroOneNorm <- function(dat){
  
  # dat <- rnorm(30,0,4)
  
  dOut <- dat - min(dat,na.rm=T)
  dOut <- dOut/max(dOut,na.rm=T)
  dOut <- dOut/(1/.998)
  dOut <- dOut + .001
  
  return(dOut)
}


