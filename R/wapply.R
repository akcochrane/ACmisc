
#' Apply a function to windows of the data
#' 
#' Iteratively applies function \code{FUN} to subsets of the data defined by a 
#' sliding window centered on each index in turn.
#' 
#' For each index, the indices \code{current_index - floor(wsize/2)} to 
#' \code{current_index - floor(wsize/2)} are taken and the function \code{FUN}
#' is applied to that subset of data. At the edges, where \code{floor(wsize/2)}
#' would extend beyond the data, the out-of-range indices are ignored. This means 
#' that at the edges of the data the window size is effectively around half of 
#' the full window size.
#' 
#' The values returned from the function are simplified into a vector if possible,
#' otherwise to a matrix, and as a last resort as a list.
#' 
#' As this method iterates over each index (in a matrix or data frame, each row),
#' it is of the utmost importance that the data is sorted appropriately (which is
#' equally applicable for similar functions such as \code{runmed()}).
#'
#' @param x Data; vector, matrix, or data frame.
#' @param wsize Window size, in number of indices. 
#' @param FUN Function to apply to each window of data. If applying to a vector \code{FUN} can be a simple function like \code{median} (see example); otherwise it should be defined as having one argument, \code{x}, indicating data.
#'
#' @export
#'
#' @examples
#' 
#' d <- rnorm(10)
#' window_median_runmed <- runmed(d,3) # from `stats` package
#' window_median_wapply <- wapply(d, 3, median) # identical except at edges
#' 
#' d <- data.frame(x = rnorm(400))
#' d$y <- d$x + rnorm(400)*(sin(seq(0,12,length = 400)+1))
#' 
#' window_correl <- wapply(d, 40, function(x){cor(x[,'x'],x[,'y'])})
#' 
#' plot(window_correl)
#' 
wapply <- function(x,wsize,FUN){
  
  wsizeHalf <- max(c(1,floor(wsize/2)))
  
  if(is.vector(d)){ndim <- 1 ; try({vectLength <- length(x)},silent=T)}
  try({if(length(dim(d))==2){ndim <- 2; vectLength <- nrow(x)}},silent=T)
  
  
  outList <- list()
  for(curCase in 1:vectLength){
    startInd <- max(c(1,curCase - wsizeHalf))
    endInd <- min(c(vectLength,curCase + wsizeHalf))
    if(ndim == 1){
      outList[[curCase]] <- 
        do.call('FUN',args = list(x = x[startInd:endInd]))
      
    }
    if(ndim == 2){
      outList[[curCase]] <- 
        do.call('FUN',args = list(x = x[startInd:endInd,]))
    }
  }
  
  outLengths <- lapply(outList,length)
  
  if(all(outLengths == 1)){ # return a vector
    dOut <- unlist(outList)
    return(dOut)
  }
  
  if(length(unique(outLengths))==1){ # return a matrix
    try({
      dOut <- do.call('cbind',outList)
      return(dOut)
    },silent=T)
  }
  
  return(dOut) # return a list
}
