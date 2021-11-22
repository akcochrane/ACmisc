#' Wilcoxon tests, returning a z value that respects the sign of the differences between two vectors
#'
#' Essentially, a wrapper for \code{ wilcox.test} that also returns the p value, 
#' Z value, W value, and the difference in medians.
#'
#' @param x First vector to compare
#' @param y Second vector to compare
#' @param pairTest Are x and y paired?
#'
#' @export
#' 
#' @examples 
#' # Note that this doesn't make *a lot* of sense, since they're two totally different measurements
#' wilcZ(dat_cochraneEtAl_2019_PLOSOne$ufov, dat_cochraneEtAl_2019_PLOSOne$mot,pairTest = T)
#' 

wilcZ <- function(x,y,pairTest=F){
  
  curTest <- wilcox.test(x,y,paired=pairTest)
  medDiff <- median(y,na.rm=T)-median(x,na.rm=T)
  
  return(data.frame(z=qnorm(curTest$p.value)*sign(medDiff),
                    p=curTest$p.value,
                    w=curTest$statistic,
                    medDiff = medDiff
                      ))
}
