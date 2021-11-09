#' Re-randomize seed
#'
#' If the seed has been set, this will un-set it.
#'
#' @export
#' 
#' @examples 
#' set.seed(13) # set an arbitrary seed
#' rnorm(5)     # see the resulting "random" values
#' set.seed(13) # re-set the same arbitrary seed
#' rnorm(5)     # see the same "random" values
#' set.seed(13) # re-set the same arbitrary seed
#' resetSeed()  # un-set the current seed
#' rnorm(5)     # see different "random" values
#'
resetSeed <- function(){

  rm(.Random.seed, envir=globalenv())
  
  # rngHole <- rnorm(as.numeric(substr(format(Sys.time(),'%s'),7,10)))
  # rm(rngHole)
  
}
