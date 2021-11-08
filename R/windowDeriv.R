
#' Get by-time Level, Change and Acceleration
#'
#' Fits a quadratic OLS regression model to each item in \code{dat}. Iteratively centers \code{tim} on each item
#' in \code{dat}, then weights the data using a Gaussian, with the half width half max (\code{hwhm}), around that \code{tim}.
#'
#' Probably the most useful part of resampling the data is to get the relative uncertainty of estimates,
#' to possibly reject the estimates with the highest uncertainty of estimates.
#'
#' To do:
#'
#' - How to effectively orthogonalize d0 and d2?
#'
#' - Is this robust to non-sorted times? Probably not
#'
#' - How to recommend appropriate HWHM? Does the current implementation work
#'
#' @param dat Vector of values to get time-dependent derivatives. Best to sort by \code{tim}.
#' @param tim Vector of time values. Should be positive real numbers. Best to be sorted.
#' @param hwhm Distance (in time) at which values should receive half of the weight as the centered time. Approximately \code{1.18*sd} of the Gaussian weights. Defaults to 2 times \code{mean(diff(tim))}
#' @param detrend Logical. Should an exponential trend, with a half-life of .5*max(tim), be subtracted prior to fitting derivatives?
#' @param orthog NOT IMPLEMENTED... always F for now. in the future, use contrasts to make derivatives orthogonal?
#' @param nResample Number of resamples, to get SE in estimates
#' @param resamplePercent Percent of data to resample, to get pseudo-SE in estimates. *this is currently arbitrary, and the scaling pseudo-SE is therefore arbitrary as well. However, relative pseudo-SE should be interpretable.*
#' @export
#'
#' @examples
#' d <- data.frame(tim = 1:100, y = rnorm(100) + log(seq(1,100))) ; d$y[sample(100,5)] <- NA
#' fit <- windowDeriv(d$y,d$tim,nResample = 20) # should probably resample >100 to get anything like reliable estimates
#' plot(d) ; lines(d$tim,fit$mean_est$level)
#' psych::pairs.panels(fit$mean_est)
#'
#' # plot relative uncertainty in estimates
#' psych::pairs.panels(data.frame(Time=1:nrow(d),fit$resampled_se_est[,2:4]))
#'
#' ## look at frequency estimates for sunspot data
#' dSun <- data.frame(spots = as.numeric(datasets::sunspot.month))
#' dSun$month <- 1:nrow(dSun) / 12
#' # takes a while to run:
#' plot(c(0,2),c(0,20),col = 'white',xlab='hwhm (years)',ylab='estimated sunspot phase length',sub='true phase length is approximately 11')
#' for (curHWHM in seq(.1,2,.1)){points(curHWHM,attr(windowDeriv(dSun$spots,dSun$month, hwhm = curHWHM),'frequency'))}
#' for (curHWHM in seq(.1,2,.1)){points(curHWHM,attr(windowDeriv(dSun$spots,dSun$month, hwhm = curHWHM,orthog=T),'frequency'),col='red')}
#'
windowDeriv <- function(dat,tim,hwhm='auto_2',detrend=F,orthog=F,nResample = 0, resamplePercent = .75){

  hwhmScale <- mean(diff(tim),na.rm=T)
  if(hwhm == 'auto_2'){
    hwhm <- hwhmScale*2
  }

  SD = hwhm * (1/sqrt(2*log(2))) # convert HWHM to SD

  nDat <- length(dat)

  ## detrend using a saturating exponential, if desired
  if(detrend){
    expoTime <- 2^((min(tim,na.rm=T)-tim)/(max(tim,na.rm=T)/2))
    trendLM <- lm(dat ~ expoTime)
    dat <- dat - (trendLM$coefficients[1] + trendLM$coefficients[2]*expoTime)
    rm(expoTime,trendLM)
  }

  timLM <-  matrix(NA,nDat,3) ; for(curIndex in 1:nDat){

    ## ## be robust to NAs in the time vector (which shouldn't really happen)
    if(!is.na(tim[curIndex])){ ## if time was not NA
      curTime <- tim[curIndex]
    }else{curTime <- mean(tim[(curIndex-1):(curIndex+1)],na.rm = T)} ## if point time is NA, then try using 3-trial average

    ## ## fit models
    if(!is.na(curTime)){## if time was not NA
      curDat <- data.frame(dat=dat,tim = tim - curTime)
      curDat$tim2 <- (curDat$tim ^ 2)/2

      modWeight <- dnorm(tim, mean=tim[curIndex], sd = SD)

      if(!orthog){ ## NOT IMPLEMENTED

        timLM[curIndex,] <- lm(dat ~ tim + tim2, curDat, w=modWeight )$coefficients

      }else{ ## DOES NOT PROVIDE ORTHOGONAL ESTIMATES

        ## level
        {
          curLM0 <- lm(dat ~ 1,  curDat, w= modWeight)
          timLM[curIndex,1] <- curLM0$coefficients
        }

        ## velocity
        {
          curDat$dat_c0 <- residuals(curLM0) # should do this by hand, to be robust to NA
          curLM1 <- lm(dat_c0 ~ 0 + tim, curDat,  w=modWeight)
          timLM[curIndex,2] <- curLM1$coefficients
        }

        ## acceleration
        {
          curDat$dat_c01 <- residuals(lm(dat ~ 1 + tim, curDat, w=modWeight)) # should do this by hand, to be robust to NA
          curLM2 <- lm(dat_c01 ~ 0 + tim2, curDat, w=modWeight)
          timLM[curIndex,3] <- curLM2$coefficients
        }
      }
    }
  }





  mean_est = data.frame(
    tim = tim,
    dat = dat,
    level = timLM[,1],
    velocity = timLM[,2],
    acceleration = timLM[,3]
  )

  mod_damped_linear_oscillator <- lm(acceleration~level+velocity,mean_est)

  attr(mean_est,'DLO') <- mod_damped_linear_oscillator
  try({ ## in case there is a pathological (positive) relationship between d0 and d2
    suppressWarnings({phase_length <- 2*pi/sqrt(-mod_damped_linear_oscillator$coefficients['level'])})
    names(phase_length) <- 'phase_length'
    attr(mean_est,'frequency') <- phase_length
  },silent=T)

  if(nResample == 0){
    return(mean_est)
  }else{

    ## ##
    ## Get resampled [subsampled] datasets to re-run the analysis
    bootLM <-  array(NA,c(nDat,5,nResample))

    for (curBoot in 1:nResample){
      tmpDat <- dat
      tmpDat[sample(length(dat),round(length(dat)*resamplePercent))] <- NA
      tmpMod <- as.matrix(windowDeriv(dat=tmpDat,tim=tim,hwhm=hwhm,orthog=orthog,nResample=0))
      bootLM[,,curBoot] <- tmpMod

    }

    resampled_se_est <- data.frame((
      apply(bootLM,1:2,function(x){quantile(x,c(.975),na.rm=T)}) -
        apply(bootLM,1:2,function(x){quantile(x,c(.025),na.rm=T)})
    )/(2*qnorm(.975)))

    colnames(resampled_se_est) <- colnames(mean_est)


    return(list(mean_est=mean_est,
                resampled_se_est = resampled_se_est,
                resampled_array = bootLM))
  }
}
