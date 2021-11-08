

#' Fit a Quick [see Ahissar & Hochstein] function, given data [currently nonsense fits]
#'
#' Given data, what are the predicted Quick parameters?
#' Fits [by minimizing negative log likelihood] either a time-evolving or a static Quick function 
#' (Ahissar & Hochstein). If time-evolving,
#' then threshold is each fit by a 3-parameter exponential decay function.
#'
#' Returns a data frame with fit parameters and the [negative] log likelihood.
#'
#' @param x stimulus strengths to fit [positive, e.g., SOA]
#' @param y responses to fit [bernoulli]
#' @param t vector of trial numbers, if desired
#' @param model "static" or "exp," depending on what the parameters correspond to
#' @param lapseRate At what lapse rate were the parameters fit?
#' @param maxThresh How high of a threshold is acceptable (in multiples of the largest stimulus strength value)?
#' @keywords  keywords
#' @export
#' @examples
#' example

fitQuickfun <- function(x, y, t=NA, model='static',lapseRate = .02,numGuesses=25,maxThresh=5){

quickfunError<- function(p, x, y, t=NA, model='static',prob=.75,minRate=0,lapseRate = .02){
  y1<- quickfun(p, x, t, model,lapseRate=lapseRate,minRate=minRate)*.998+.001
  L = sum(y*log(y1) + (1-y)*log(1-y1))
  error<- -L
  if(max(abs(y1-.5),na.rm=T)<.1){error=error^2} # this is to make the error spike, if the model never predicts a probability of greater than .6
  return(error)
}

# #

if (model=='static'){
  guesses <- cbind(runif(numGuesses,0,max(abs(x),na.rm=T)), # threshold guesses
                   runif(numGuesses,0,max(abs(x),na.rm=T))# slope guesses
  )
  fits <- data.frame(th=rep(NA,numGuesses),
                     slo=NA,
                     err=NA)
}

if (model=='exp'){
  guesses <- cbind(runif(numGuesses,0,max(abs(x),na.rm=T)), # threshold asymptote guesses
                   runif(numGuesses,0,max(abs(x),na.rm=T)), # threshold scale guesses
                   runif(numGuesses,1,log2(length(x))+1), # threshold rate guesses

                   runif(numGuesses,0,max(x,na.rm=T)) # slope guesses
                   )
  fits <- data.frame(thAs=rep(NA,numGuesses),
                   thSt=NA,
                   thRa=NA,
                   slope=NA,
                   err=NA)
}



for (fitNum in 1:numGuesses){
  try({
curFit<- optim(guesses[fitNum,],
               fn=quickfunError,
               x=x,
               y=y,
               t=t,
               model=model,
               lapseRate=lapseRate,
               minRate=2,
               lower=c(0,0,1,0),
               upper=c(max(x,na.rm=T)*maxThresh,max(x,na.rm=T)*maxThresh,log2(length(x))+5,3),
               method="L-BFGS-B",
               control=list(maxit=10000))
cat('. ')

fits[fitNum,] <- c(curFit$par,curFit$value)
rm(curFit)
},silent=F)
}
fits <- na.omit(fits)

bestFit <- fits[which.min(fits$err)[1],]

# predVals <- quickfun(unlist(bestFit)[1:4],x=x,t=t) ; runningAvg(predVals,toPlot = T)

return(bestFit)

}
