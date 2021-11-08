

#' Put lmer results into a data frame
#'
#' Input an lmer object, and up to three optional descriptions. In return, you get a table
#' with the model's estimates of F, p, df, b, and CIs.
#'
#' Unfortunately, this also creates a modelOutputByproducts.txt file... Which is a compromise,
#' to prevent lmSupport from blowing up the console with a ton of output.
#'
#' @param modIn lmer object
#' @param modID1 optional identification variable (e.g., dataframe name)
#' @param modID2 optional identification variable (e.g., task name)
#' @param modID3 optional identification variable (e.g., a grouping variable)
#' @keywords  keywords
#' @export

lmer2table<- function(modIn,modID1=NA,modID2=NA,modID3=NA,round2=4){

  require(lmSupport)
  options(warn=-1)
  sink('modelOutputByproducts.txt')
  mSumm <- modelSummary(modIn,t=F,Print=F)
  sink(NULL)
  options(warn=0)

  outTable <- data.frame(
    modID1=modID1,
    modID2=modID2,
    modID3=modID3,
    DV=substr(as.character(mSumm$call)[2],1,4),
    IV=rownames(mSumm$KRAppox),
    fval=mSumm$KRAppox[,'F'],
    pval=mSumm$KRAppox[,'Pr(>F)'],
    df=mSumm$KRAppox[,'error df'],
    estim=mSumm$KRAppox[,'Estimate'],
    #
    # CIs using T distribution given the KR approximated df:
    ciLo=mSumm$KRAppox[,'Estimate']-mSumm$KRAppox[,'SE']*qt(.975,mSumm$KRAppox[,'error df']),
    ciHi=mSumm$KRAppox[,'Estimate']+mSumm$KRAppox[,'SE']*qt(.975,mSumm$KRAppox[,'error df']),
    fullCall=as.character(modIn@call)[2]
  )

  if(is.numeric(round2)){
    outTable[,6:11] <- round(outTable[,6:11],round2)
  }

  return(outTable)
}

if(F){
  library(AaronFuns) ; library(lme4)

  dat <- feedbackRawDat[1:2000,]
  my_mod <- lmer(resp~stimStrength+(stimStrength|subid),data=dat)
  lmer2table(m)
}
