

#' Find a robust bootstrapped correlation
#'
#' First optimally transforms x and y using \code{\link{YeoJohn}} (i.e., Yeo-Johnsons transformation
#' with a lambda optimized to minimize skew).
#' Then resamples the data 2000 times and, for each resample, calculates the
#' Pearson correlation coefficient along with the corresponding T value and
#' log-3 Bayes Factor (calculated with a default \code{BayesFactor::correlationBF(x,y,'medium')}).
#' Returns a formatted string \code{$results},
#' the full resampled set of values in \code{$resampled_tests}, as well as
#' individual values for the median \code{$r} and its corresponding \code{$t} and log-2 \code{$bf}.
#' The \code{$data} is returned in a data frame with variables \code{xt} [transformed x] and
#' \code{yt} [transformed y], each of which has a set of attributes (e.g., the "rawDat" and the Yeo-Johnson "lambda").
#'
#' The 'medium' prior for \code{BayesFactor::correlationBF} is a transformed beta(3,3) (
#' see the BayesFactor documentation as well as Ly, Verhagen, and Wagenmakers (2015) )
#'
#' Note that this function will not work with fewer than 9 cases. Because you probably shouldn't be
#' calculating correlations with only 8 cases.
#'
#' @param x numeric vector
#' @param y numeric vector
#'
#' @export
#'
robustCor <- function(x,y){

  xt <- YeoJohn(x)
  yt <- YeoJohn(y)

  dft <- data.frame(xt = xt, yt = yt)

  allCorr <- replicate(2000,{

    corDF <- 2 ; while(corDF < 9){ # prevent pathological resampling
    curSample <- na.omit(dft[sample(length(x),replace = T),])
    corDF <- nrow(curSample) - 2
    }

    corR <- cor(curSample$xt,curSample$yt,method='pearson')

    corT <- psych::r2t(corR,corDF)
    corBF <- suppressMessages({log(exp(
      BayesFactor::correlationBF(curSample$xt,curSample$yt)@bayesFactor$bf
      ),base=3)})
    names(corBF) <- 'log3(BF10)'

    c(corR = corR, corDF = corDF, corT = corT, corBF = corBF)

  })

  corR <- median(allCorr['corR',])
  corT <- allCorr['corT',which.min(abs(allCorr['corR',] - corR))]
  corBF <- allCorr['corBF.log3(BF10)',which.min(abs(allCorr['corR',] - corR))]
  corDF <- allCorr['corDF',which.min(abs(allCorr['corR',] - corR))]


    return(list(results = paste0('r(',corDF,') = ',
                                 round(corR,2),
                                 ' [',
                                 paste(round(quantile(allCorr['corR',],c(.025,.975)),2),collapse=','),
                                 '], BF~log3~ = ',
                                 round(corBF,2))
                ,r = corR
                ,t = corT
                ,bf = corBF
                ,resampled_tests = data.frame(t(allCorr))
                ,data = dft
    ))
}
