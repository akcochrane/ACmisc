#' Estimate matrix of bivariate out-of-sample variance explained
#'
#' Uses \code{tef_rlm_boot()} to get out-of-sample variance explained as well as bootstrapped p values.
#' The latter is two times the proportion of the bootstrapped distribution is on the opposite side of
#' zero than the median coefficient.
#'
#'
#' @param xvars Vector of variable names, with which to estimate bivariate relations.
#' @param partialvars If a vector of variable names are provided, then the lower triangle of the output is the out-of-sample variance explained \emph{when controlling for the variables in \code{partialvars}}.
#' @param data Data to use
#' @param nBoots Number of resamples
#'
#' @export
#'
oos_partial <- function(xvars,partialvars = c(),data,nBoots = 200){
  library(MASS)
  library(TEfits)
  library(ggplot2)

  if(F){ ## for testing
    xvars <- c('x1','x2','x3')
    partialvars <- c('p1','p2')
    data <- data.frame(x1 = rnorm(100))
    data$x2 <- data$x1 + rnorm(100)
    data$x3 <- data$x1+ rnorm(100,0,.5)
    data$p1 <- data$x2 + rnorm(100,0,.01)
    data$p2 <- data$x3 + rnorm(100,0,.01)
  }

  dOut <- matrix(NA
                 ,nrow  = length(xvars)
                 ,ncol = length(xvars)
                 ,dimnames = list(xvars,xvars)
  )
  dOut_bootP <- matrix(NA
                       ,nrow  = length(xvars)
                       ,ncol = length(xvars)
                       ,dimnames = list(xvars,xvars)
  )

  for(curVar1 in 1:length(xvars)){
    for(curVar2 in 1:length(xvars)){

      ## ## Upper triangle is raw bivariates
      if(curVar1 < curVar2){
        rlm1 <- tef_rlm_boot(
          formula(paste(xvars[curVar1],'~',xvars[curVar2]))
          ,datIn = data, nBoot = nBoots
        )
        rlm2 <- tef_rlm_boot(
          formula(paste(xvars[curVar2],'~',xvars[curVar1]))
          ,datIn = data, nBoot = nBoots
        )
        dOut[curVar1,curVar2] <-
          mean(c(
            rlm1$bootSummary[xvars[curVar2],'dRsq_oos']
            ,rlm2$bootSummary[xvars[curVar1],'dRsq_oos']
          ))
        dOut_bootP[curVar1,curVar2] <-
          mean(c(
            rlm1$bootSummary[xvars[curVar2],'bootP']
            ,rlm2$bootSummary[xvars[curVar1],'bootP']
          ))
        suppressWarnings(rm(rlm1,rlm2))
      }

      ## ## Lower triangle is partial bivariates
      if(curVar1 > curVar2 && length(partialvars)>0){
        rlm1 <- tef_rlm_boot(
          formula(paste(xvars[curVar1],'~',xvars[curVar2],'+',
                        paste(partialvars,collapse = '+')))
          ,datIn = data, nBoot = nBoots
        )
        rlm2 <- tef_rlm_boot(
          formula(paste(xvars[curVar2],'~',xvars[curVar1],'+',
                        paste(partialvars,collapse = '+') ))
          ,datIn = data, nBoot = nBoots
        )
        dOut[curVar1,curVar2] <-
          mean(c(
            rlm1$bootSummary[xvars[curVar2],'dRsq_oos']
            ,rlm2$bootSummary[xvars[curVar1],'dRsq_oos']
          ))
        dOut_bootP[curVar1,curVar2] <-
          mean(c(
            rlm1$bootSummary[xvars[curVar2],'bootP']
            ,rlm2$bootSummary[xvars[curVar1],'bootP']
          ))
        suppressWarnings(rm(rlm1,rlm2))
      }


      if(curVar1 == curVar2){
        # nothing yet
      }

    }
  }

  return(list(
    dRrsq_oos = dOut
    ,bootP = dOut_bootP
  ))
}
