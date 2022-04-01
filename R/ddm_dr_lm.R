#' Fit linear predictors of DDM drift rate
#'
#' Essentially fits a generalized linear model for a Wiener process, using the density
#' estimation function from the \code{RWiener} package. First fits a standard 4-parameter Wiener
#' model to the entire RT+response-category vector. Then uses the boundary separation, bias, and
#' non-decision time from this overall model, and finds the best set of parameters to create a
#' drift rate vector as a linear function of the right-hand side of \code{formIn}. \strong{Only
#' numeric predictors have been tested}.
#' 
#' Note that, for parameters other than drift rate, the \code{RWiener} labeling conventions are used: 
#' \code{alpha} is boundary separation, \code{tau} is non-decision time, and
#' \code{beta} is bias.
#' 
#' \code{rt ~ 1} will cause an error; in principle, if a unconditional point estimate of drift rate is desired,
#' a dummy constant numeric variable (e.g., all 0s) could be included as a predictor.
#'
#' @param formIn formula. Left hand side should be response time.
#' @param data data
#' @param respVar name of the "boundary" column (argument is character; data must be binary 0 and 1)
#' @param nBoots number of bootstrapped samples for CI of drift rate parameters. Must be an integer >1
#' @param fixBias Optional scalar value to fix the bias [beta] parameter of the Wiener function
#'
#' @seealso 
#' \code{\link[RWiener]{wdm}} for the method upon which this is based; 
#' \code{\link[brms]{brm}} with family "wiener" for a much better (but slower) option.
#'
#' @return
#' @export
#'
#' @examples
#' d <- data.frame(rt = .2 + exp(rnorm(50,-.5,.5))
#' ,corr = rbinom(50,1,.5)
#' ,totalTrialNum = 1:50
#' )
#' d$stimStr <- d$rt+rnorm(50,d$corr)
#' m <- ddm_dr_lm(rt ~ totalTrialNum * stimStr, data = d,respVar='corr',nBoots = 30)
#' m$model
#'
ddm_dr_lm <- function(formIn,data,respVar,nBoots=2,fixBias = NA){
  require(RWiener)
  
  # to do: option between fitting ndt etc for each bootstrap, or keeping the same one for overall

  yVar <- data[,as.character(formIn[[2]])] * (2*data[,respVar] - 1)

  failed = T ; nTries <- 0 ; while(failed && nTries < 100){
    nTries <- nTries+1
    try({
  if(!is.numeric(fixBias)){
  total_wdm <- wdm(yVar)}else{
  total_wdm <- wdm(yVar,beta = fixBias)}
  failed <- F},silent=T)
  }
  if(failed){stop('RWiener::wdm() fit failed. Please check your data and respVar.')}

  modTerms <- c('dr_Intercept')
  modRHS <- modTerms

  for(curTerm in attr(terms(formIn),"term.labels")){
    modTerms <- c(modTerms,paste0('dr_',gsub(':','_',curTerm)))
    modRHS <- paste0(modRHS,' + ',modTerms[length(modTerms)],'*',gsub(':','*',curTerm))
  }

  findDRerr <- function(pars,lhs=yVar,rhs=modRHS,modTerms=modTerms,data,alpha,tau,beta){

    tmpDat <- data.frame(t(pars),data)
    colnames(tmpDat)[1:length(pars)] <- modTerms

    dr <- eval(expr=as.formula(paste('~',modRHS))[[2]],env=tmpDat)


    dr_dens <- c()
    for (curDR in 1:length(dr)){ # for the life of me, I can't seem to vectorize dwiener()
    tryCatch({
    dr_dens <- c(dr_dens,dwiener(lhs[curDR,'q'],alpha=alpha,beta=beta,tau=tau,delta=dr[curDR],resp=lhs[curDR,'resp']))
    },
    error = function(x){ dr_dens <- c(dr_dens, 1E-20)}
    )
    }

    negLL <- -sum(log(dr_dens))
    if(is.infinite(negLL)){
      negLL <- 1E15
    }
    return(negLL)
  }


  startPars <- c(total_wdm$coefficients['delta'],rep(0,length(modTerms)-1))
  names(startPars)<- modTerms

  fitMod_all <- optim(startPars,
        findDRerr,
        lhs=revamp(as.wiener(yVar)),
        rhs=modRHS,
        modTerms=modTerms,
        data=data,
        alpha=total_wdm$coefficients['alpha'],
        tau=total_wdm$coefficients['tau'],
        beta=total_wdm$coefficients['beta'])


  fitMods_boots <- as.data.frame(t(replicate(nBoots,
                             {
                               fitPars <- NA
                               while(!is.numeric(fitPars)){ # protects from pathological sampling
                               curDat <- sample(nrow(data),replace=T)
                               fitPars <- optim(startPars,
                                     findDRerr,
                                     lhs=revamp(as.wiener(yVar[curDat])),
                                     rhs=modRHS,
                                     modTerms=modTerms,
                                     data=data[curDat,],
                                     alpha=total_wdm$coefficients['alpha'],
                                     tau=total_wdm$coefficients['tau'],
                                     beta=total_wdm$coefficients['beta'])$par
                               }
                               names(fitPars) <- modTerms
                               fitPars
                             })))

  model <- data.frame(Estimate=c(fitMod_all$par,
                                         total_wdm$coefficients['alpha'],
                                         total_wdm$coefficients['tau'],
                                         total_wdm$coefficients['beta']),
                                         ci025 = NA , ci975 = NA , nBoot = NA)

  attr(model,'neg_LL') <- fitMod_all$value
  
  model[1:length(modTerms),'ci025'] <- apply(fitMods_boots,2,quantile,.025)
  model[1:length(modTerms),'ci975'] <- apply(fitMods_boots,2,quantile,.975)
  
  
  model[1:length(modTerms),'nBoot'] <- nBoots
  
  tmpDat <- data.frame(t(fitMod_all$par),data)
  colnames(tmpDat)[1:length(fitMod_all$par)] <- modTerms
  
  data$fitted_DR <- eval(expr=as.formula(paste('~',modRHS))[[2]],env=tmpDat)
# 
#   fit_DR <- rowSums(
#     sweep(m_ddm$data[,basisVars] , 
#           2,
#           m_ddm$model[paste0('dr_',basisVars),'Estimate'],
#           '*'
#     )
#   )
  
  return(list(model=model,wdmFit = total_wdm,pars_boots = fitMods_boots,data=data))
}
