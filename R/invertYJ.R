
#' Invert a Yeo-Johnson transformation
#'
#' @param xt Transformed values of data
#' @param lambda Lambda value that had transformed data into xt
#' @param interval Interval in which to look for data. If search is difficult, restricting the search to plausible values should help.
#' @param keepMedian Logical. Should the original data's median be kept?
#' @param keepMAD Logical. Should the original data's MAD be kept?
#'
#' @export
#'
#' @examples
#' x <- 1:10
#' x1 <- VGAM::yeo.johnson(x,2)
#' x_recovered <- invertYJ(x1,2)
#'
#' d <- rnorm(5,3) + exp(5)
#' d1 <- AaronFuns::YeoJohn(d)
#' invertYJ(attr(d1,'transformedDat'),attr(d1,'lambda'),interval = c(-1e3,1E3)) # fails with too large of an interval.
invertYJ <- function(xt, lambda, interval = c(-1E10,1E10), keepMedian = T, keepMAD = T){

invertYJerror <- function(x,xt, lambda){sum((xt - VGAM::yeo.johnson(x,lambda))^2)}

x <- c() ;for(curXT in xt){
tryCatch({
  x <- c(x,round(optimize(invertYJerror,interval=interval,xt=curXT,lambda = lambda,tol = 1E-7)$minimum,6))}
,error = function(.){x <- c(x,NA)}
)

}
  return(x)
}
