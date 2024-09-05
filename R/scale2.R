
#' Scale a variable to a new mean and standard deviation
#'
#' @param x Data to scale
#' @param newmean Desired mean
#' @param newsd Desired standard deviation
#' @param yj [logical] Whether Yeo-Johnson transformation should be applied to minimize univariate skew
#'
#' @export
#'
#' @examples
#' 
#' y <- c(rnorm(200,3,3),rnorm(100,10,3) )
#' hist(y, breaks = 40)
#' hist(scale2(y,4,2), breaks = 40)
#' hist(scale2(y,-1,.5), breaks = 40)
#' hist(scale2(y,10,1,yj=T),breaks = 40)
#' hist(scale2(exp(y),10,1,yj=F),breaks = 40)
#' hist(scale2(exp(y),10,1,yj=T),breaks = 40)
#' 
scale2 <- function(x,newmean = 0, newsd = 1, yj = F){
  
  if(yj){x <- ACmisc::YeoJohn(x)}
  
  x_mean <- mean(x, na.rm=T)
  x_sd <- sd(x, na.rm = T)
  x_scaled <- (x - x_mean)/x_sd
  
  return(x_scaled * newsd + newmean)
  
}