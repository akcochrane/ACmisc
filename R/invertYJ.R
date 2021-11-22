
#' Invert a Yeo-Johnson transformation
#' 
#' Once a Yeo-Johnson transformation has been applied, is may be useful to un-transform the data
#' to recover the raw data. See examples.
#' 
#' Because \code{\link{YeoJohn}} defaults to retaining the default median and MAD of the original data,
#' it is not trivial to simply apply the Yeo-Johnson lambda "in reverse." Instead, this function
#' uses 1-dimensional OLS optimization to find the vector of data that has the same median and mean as 
#' the transformed data, but also has the lambda applied "in reverse." Such optimization has limits, however;
#' both \code{\link{YeoJohn}} and \code{\link{invertYJ}} will work best on relatively small values (e.g., absolute 
#' values less than 1000) and on datasets with more than 10 numbers or so.
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
#' d <- rnorm(10,3) + rexp(10,4)
#' d1 <- ACmisc::YeoJohn(d)
#' invertYJ(attr(d1,'transformedDat'),attr(d1,'lambda'),interval = c(-1e3,1E3)) # fails if there is too large of an interval.
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
