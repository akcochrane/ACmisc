#' Wilcoxon tests, returning a z value that respects the sign of the differences between two vectors
#'
#' Also returns the p value, W value, and the difference in medians.
#'
#' @param x First vector to compare
#' @param y Second vector to compare
#' @param pairTest Are x and y paired?
#'
#' @export

wilcZ <- function(x,y,pairTest=F){
  
  curTest <- wilcox.test(x,y,paired=pairTest)
  medDiff <- median(y,na.rm=T)-median(x,na.rm=T)
  
  return(data.frame(z=qnorm(curTest$p.value)*sign(medDiff),
                    p=curTest$p.value,
                    w=curTest$statistic,
                    medDiff = medDiff
                      ))
}
