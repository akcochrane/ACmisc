

#' Find Lower RT bound, return the bound and a vector of "keep values"
#'
#' Following the recommendation of Ratcliff & Tuerlinckx (2002)
#' and finds the threshold at which, below that RT, responses are not different than chance
#' (binomial test p-value > .05). RT thresholds are tested in steps of .001 if RT is 
#' in seconds, steps of 1 if RT is in milliseconds. Begins by rejecting the lowest
#' 10 RTs to start the binomial tests. If the number rejected is below 10 or so, 
#' consider not rejecting any! Returns lowerRT,
#' the rejection level, as well as keepVect, a vector tagging which RTs to keep. 
#'
#' @param rtVector vector of RTs
#' @param accVector vector of 0 and 1 or T and F indicating correct performance
#' 
#' @export
#' 
findLowerRT <- function(rtVector,accVector){
  
  # NOTE: 
    # SHOULD EDIT FUNCTION TO START AT 10 RTS AND IF THIS NUMBER IS GOOD, RETURN 
    # A REJECTION THRESHOLD OF 0 AND A VECTOR OF ALL ACCEPTANCE.
  
  # start by rejecting the trials with the lowest 4 RTs
  ## ## if the number rejected is below 10 or so, consider not rejecting any!
  rejectLevel  <- sort(rtVector)[10] 
  
  ## if RT is in ms, step should be 1. if it is in seconds, step should be .001
  stepSize <- 10^floor(log10(rejectLevel)-2)
  if(stepSize<.001){stepSize <- .001}
  
  reject <- F
  while (reject==F){
    
    relevantData <- accVector[rtVector < rejectLevel]
    
    if(binom.test(sum(relevantData),length(relevantData),p=.5,alternative='greater')$`p.value` < .05 ){
      reject=T
    }else{rejectLevel <- rejectLevel + stepSize}
  }
  return(list(lowerRT=rejectLevel,keepVect = rtVector > rejectLevel))
}
