
#' Maximum-likelihood fitting of a brm-style model
#'
#' @param formula Model formula
#' @param data Model data
#' @param ... Additional arguments to pass to an initial \code{brm} model
#' @param iter Number of iterations. If \code{strata == NA} then this is draws from a normal approximation (inverse Hessian); if \code{strata != NA} then this is the number of bootstrap samples. 
#' @param brmModel An initial brm model to re-fit using maximum-likelihood.
#' @param strata The strata to use in resampling; a vector the length of the data set's rows. If all values are equal, no stratification is used in resampling
#' @param crossvals Number of cross-validations [currently not implemented]
#' @param bootstrap Bootstrapping is done by default if \code{strata} is provided. However, in certain cases (e.g., if cross-validation is desired but bootstrapping is not), bootstrapping can be turned off by setting \code{bootstrap = F}
#' @param refresh Should progress be printed?
#' 
#' Use \code{\link[rstan]{optimizing}} fit a \code{\link[brms]{brm}}-style model.
#' 
#' This function is an attempt to merge the usability and flexibility of \code{\link[brms]{brm}}
#' and the [potential] efficiency of \code{\link[rstan]{optimizing}}. This may be useful in various
#' circumstances, for example, when desiring a direct comparison between Bayesian and maximum-likelihood
#' methods (e.g., to give an indication of the influence of priors), or when other methods are
#' problematically computationally expensive. While it is always recommended to use 
#' \code{\link[brms]{brm}} with \code{algorithm = 'sampling'}, preferably checking robustness to 
#' priors and model specification through model comparison, sometimes this ideal is not going to 
#' be feasible. The current function provides maximum-likelihood estimation with a
#' \code{\link[brms]{brm}}-like interface.
#' 
#' Several major limitation persists, however. One limitation is that sampling methods 
#' inherently provide distributions of parameters whereas maximum-likelihood methods do not.
#' The simplest solution to this would be using the Hessian of the optimization, but this matrix
#' is not always well-behaved. A potentially more robust solution is to bootstrap the model (i.e.,
#' in principle, draw parameters from their sampling distribution by resampling data with replacement
#' and re-fitting the model), but this leads to the second limitation. 
#' As of now the implementation of \code{brm_optimizing}
#' is apt to crash the R session with repeated calls or with very large models/datasets; 
#' bootstrapping as well as robustness checks of point estimates are necessary to fully
#' interpret maximum-likelihood fits, but these are limited when R is apt to crash.
#' 
#' In summary, the user is advised to proceed with caution. Repeated maximum-likelihood fitting
#' and bootstrapping are likely to be the most beneficial applications of this method, but the 
#' current implementation is somewhat unstable for these applications.
#' 
#' @note 
#' 
#' Most errors and warnings can be ignored; they are likely do due hiccups in the
#' pipeline that do not influence the final outcome. The error 
#' \code{Error in chol.default(-H)} indicates that the Hessian was ill-behaved and
#' direct approximations of the estimates' standard errors were likely impossible.
#' As long as the resulting object contains \code{Est.Error} and CI, however,
#' it is likely that at least one estimation run allowed for these approximations.
#'
#' @export
#'
#' @examples
#' 
#' m1 <- brm_optimizing(Sepal.Width ~ Sepal.Length * Petal.Width + (1|Species),iris)
#' 
#'     m_cv_1 <- brm_optimizing(Sepal.Width ~ Sepal.Length * Petal.Width + (1|Species),iris
#' ,strata = iris$Species,iter = 500, 
#' bootstrap = F, crossvals = 50)
#' m_cv_2 <- brm_optimizing(Sepal.Width ~ Sepal.Length * Petal.Width,iris
#'                          ,strata = iris$Species,iter = 500, 
#'                          bootstrap = F, crossvals = 50)
#' m_cv_1$logLik_oos - m_cv_2$logLik_oos # m_cv_1 is much better! test 1, 50 crossvals: 57.7; test 2: ; test 3: 
#' 
brm_optimizing <- function(formula=NA,data=NA,...
                           ,iter = 2000
                           ,brmModel = NA
                           , strata=NA
                           , crossvals = 0
                           , bootstrap = NA
                           , refresh = 0){
  library(brms) ; library(rstan)
  
  # to check: build boot samples initially using just row nubmers. Then save these row numbers in a matrix, so the resamples can be reconstructed
  # # but how to do this, with the bootstraps being made somewhat ad hoc piecemeal using the within-strata indices?
  
  # to test: with distributional models, make sure output names are correct
  
  # to test: all combinations of iter, strata, boostrap, and crossvals all give the desired outcomes
  
  # to check: why R crashes so often. How to avoid? Insterting gc() certainly hasn't seemed to help. 
  # # Maybe there's a Stan / Rcpp garbage collection?
  
  # to check: look up the uses and limitations of importance resampling... is it just to get the log-posterior? Or something else? 
  # Does it cause more errors, or just the same type with non-positive-definite-ness?
  
  # to check: speculate on what the sources of the differences between model fit methods are. To the extent that the a reasonable range of fit values are
  # - - - - coming out of the ML boots (e.g. sigma between .2 and 5 times the Bayesian sigma, and the distribution being gaussian-ish),
  # - - - - we can't attribute differences to pure ML instability. Especially since local minima are unlikely to be the same across boots.
  # - - - - So, the conclusion is... that HMC and priors leads to different results? possibly. What would the implications be? and WHEN?
  
  if(F){ # for testing various things
    
    m1 <- brm_optimizing(Sepal.Width ~ Sepal.Length * Petal.Width + (1|Species),iris)
    m2 <- brm_optimizing(brmsformula(Sepal.Width ~ Sepal.Length * Petal.Width + (1|Species)
                                     ,sigma ~ Petal.Width)
                         ,iris)
    combine_models(m1,m2)
    
    mb2 <- brm_optimizing(Sepal.Width ~ Sepal.Length * Petal.Width + (1|Species),iris
                          ,strata = iris$Species,iter = 20, refresh = 1)
    
    mq1 <- brm_optimizing(brmsformula(Sepal.Width ~ Sepal.Length * Petal.Width + (1|Species)) + 
                            lf(quantile = .5)
                          ,iris,family = asym_laplace(link_quantile = 'identity')) 
    
    m_cv_1 <- brm_optimizing(Sepal.Width ~ Sepal.Length * Petal.Width + (1|Species),iris
                             ,strata = iris$Species,iter = 500, 
                             bootstrap = F, crossvals = 50)
    m_cv_2 <- brm_optimizing(Sepal.Width ~ Sepal.Length * Petal.Width,iris
                             ,strata = iris$Species,iter = 500, 
                             bootstrap = F, crossvals = 50)
    m_cv_1$logLik_oos - m_cv_2$logLik_oos # m_cv_1 is much better! test 1, 50 crossvals: 57.7; test 2: ; test 3: 
    m_cv_1
    
  }
  
  # try({source('match_brm_stan_optim.R')},silent = T)
  # try({source('g:/My Drive/Aaron/functions/educational materials/match_brm_stan_optim.R')},silent = T)
  
  match_brm_stan_optim <- function(m_o_out=m_o_out,m_stan=m_stan,m_brm=m_brm,
                                   m_o_overall = m_o_overall,
                                   nSamples = nSamples
                                   ,strata=strata
                                   ,useBoot = F
  ){
    
    if(all(colnames(m_o_out) == names(m_stan))){ ## neat!
      
      for(curSampleVect in 1:length(m_stan@sim$samples[[1]])){ 
        
        m_stan@sim$samples[[1]][curSampleVect][[1]] <- 
          as.numeric(m_o_out[,names(m_stan@sim$samples[[1]])[curSampleVect]])
        
        if(names(m_stan@sim$samples[[1]])[curSampleVect] != 'lp__'){
          attr(m_stan@sim$samples[[1]],"mean_pars")[curSampleVect] <- 
            as.numeric(m_o_overall$par[names(m_stan@sim$samples[[1]])[curSampleVect]])
        }else{
          attr(m_stan@sim$samples[[1]],"mean_lp__") <- 
            as.numeric(m_o_overall$value)
        }
        
      }
      
      attr(m_stan,"args")$iter <- nSamples
      m_stan@sim$iter <- nSamples
      m_stan@sim$n_save <- nSamples
      m_stan@stan_args[[1]]$method <- 'optimizing'
      
      if(useBoot){
        if(length(unique(strata)) == 1){
          m_stan@stan_args[[1]]$algorithm <- "LBFGS + bootstrap"
        }else{
          m_stan@stan_args[[1]]$algorithm <- "LBFGS + stratified bootstrap"
        }
      }else{
        m_stan@stan_args[[1]]$algorithm <- "LBFGS + normal approximation"
      }
      
      # m_brm$m_cv <- m_cv_all
      
      m_brm$fit <- m_stan
      m_brm$strata <- strata
      
      m_brm <- rename_pars(m_brm)
      
      m_brm$max_lik_samples <- m_o_out
      
      m_brm <- add_criterion(m_brm,'bayes_R2')
      
    }else{
      stop('something went wrong with matching the ML fit to the brm output')
    }
    return(m_brm)
  }
  
  getInits <- function(m){
    lapply(m$fit@inits[[1]]
           , function(x){x + rnorm(length(x),0,.00001)})
  }
  
  useBoot <- !all(suppressWarnings(is.na(strata)))
  if(!is.na(bootstrap)){useBoot <- bootstrap}
  
  nSamples <- iter ; rm(iter)
  
  if(suppressWarnings(!is.brmsfit(brmModel))){ 
    formIn <- formula ; rm(formula)
    dataIn <- data ; rm(data)
    
    suppressMessages({suppressWarnings({
      
      m_brm <- brm(formIn,dataIn,iter = 1,chains=1
                   ,save_pars = save_pars(all = T) # is the save_pars(all = T) necessary? It might contribute to crashing for big / complex models
                   ,refresh = refresh
                   ,silent = 2
                   ,...
      )
      
    })})
    
  }else{
    suppressMessages({suppressWarnings({
      formIn <- brmModel$formula ; rm(formula)
      dataIn <- brmModel$data ; rm(data)
      m_brm <- brmModel ; rm(brmModel)
    })})
  }
  
  m_stan <- stan(model_code = m_brm$fit@stanmodel@model_code ,
                 data = make_standata(m_brm$formula ,m_brm$data)
                 ,iter=1,chains=1
                 # ,fit = m_brm$fit # would be faster, but may overflow memory
                 ,verbose = F
                 ,refresh = refresh
                 
  ) 
  
  if(length(m_stan@sim)==0){ # in case initialization failed
    
    m_stan <- stan(model_code = m_brm$fit@stanmodel@model_code
                   ,data = make_standata(m_brm$formula ,m_brm$data)
                   ,iter=1,chains=1
                   ,init_r = .01
                   ,verbose = F
                   ,refresh = refresh
                   ,init = list(getInits(m_brm))
    )
  }
  
  m_o_overall_1 <- m_o_overall_2 <- list(value = -1E20)
  nRuns <- 0; while(length(m_o_overall_1) < (1-useBoot)+4 && nRuns < 50){ try({ # the 1-useBoot thing is because the `optimizing` output with draws is 8 long and the one without draws is 4 long
    gc(verbose = F) ; Sys.sleep(.1) ; gc(verbose = F)
    m_o_overall_1 <- optimizing(m_brm$fit@stanmodel
                                ,data = make_standata(m_brm$formula ,m_brm$data)
                                # ,hessian = T # can be used for estimating the normal CI, but "draws" does it more automatically
                                ,draws = nSamples * as.numeric(!useBoot) # sets to 0 if bootstrapping is desired
                                ,iter = 10000 + nRuns*100 # increases number of iterations with each failure
                                ,as_vector=T 
                                ,init = getInits(m_brm)
    )
  },silent=(refresh==0)) ; nRuns <- nRuns + 1} ; rm(nRuns)
  #
  nRuns <- 0; while(length(m_o_overall_2) < (1-useBoot)+4 && nRuns < 50){ try({
    gc(verbose = F) ; Sys.sleep(.1) ; gc(verbose = F)
    m_o_overall_2 <- optimizing(m_brm$fit@stanmodel
                                ,data = make_standata(m_brm$formula ,m_brm$data)
                                ,draws = nSamples * as.numeric(!useBoot) # sets to 0 if bootstrapping is desired
                                ,iter = 10000 + nRuns*100
                                ,as_vector=T )
  },silent=(refresh==0)) ; nRuns <- nRuns + 1} ; rm(nRuns)
  
  if(m_o_overall_1$value > m_o_overall_2$value){
    m_o_overall <- m_o_overall_1 }else{ m_o_overall <- m_o_overall_2 }
  if(length(m_o_overall) == 1 ){stop('Model fitting failed.')}
  m_o_out <- as.data.frame(m_o_overall$theta_tilde)
  m_o_out$lp__ <-  m_o_overall$value
  
  gc(verbose = F) ; Sys.sleep(.1) ; gc(verbose = F)
  
  if(useBoot){
    nBoot <- nSamples
    smallSampleWarning <- F
    
    boot_sample_indices <- matrix(NA,nrow(m_brm$data),0)
    
    for(. in 1:nBoot){
      
      bootInds <- c() ; for(curStrata in unique(strata)){
        
        strataInds <- which(strata == curStrata)
        
        strataInds_resampled <- sample(strataInds,replace = T)
        
        bootInds <- c(bootInds,strataInds_resampled)
        
        if(length(strataInds_resampled) < 11){smallSampleWarning <- T}
        rm(strataInds,strataInds_resampled)
      }
      bootDat <- m_brm$data[bootInds,]
      boot_sample_indices <- cbind(boot_sample_indices,bootInds)
      
      m_boot_1 <- m_boot_2 <- list(value = -1E20)
      nRuns <- 0; while(length(m_boot_1) < 4 && nRuns < 50){ try({
        m_boot_1 <- optimizing(m_brm$fit@stanmodel
                               ,data = make_standata(m_brm$formula ,bootDat)
                               ,iter = 5000 + nRuns*100
                               ,as_vector=T 
                               ,init = getInits(m_brm) # for one of the two, use the brm inits
        )
      },silent=T) ; nRuns <- nRuns + 1} ; rm(nRuns)
      #
      nRuns <- 0; while(length(m_boot_2) < 4 && nRuns < 50){ try({
        m_boot_2 <- optimizing(m_brm$fit@stanmodel
                               ,data = make_standata(m_brm$formula ,bootDat)
                               ,iter = 5000 + nRuns*100
                               ,as_vector=T )
      },silent=T) ; nRuns <- nRuns + 1} ; rm(nRuns)
      
      if(m_boot_1$value > m_boot_2$value){
        m_boot <- m_boot_1 }else{ m_boot <- m_boot_2 }
      
      m_boot_out <- as.data.frame(m_boot$theta_tilde)
      m_boot_out$lp__ <-  m_boot$value
      
      m_o_out <- rbind(m_o_out,m_boot_out)
      
      rm(m_boot_out,m_boot,m_boot_1,m_boot_2)
      gc(verbose = F) ; Sys.sleep(.1) ; gc(verbose = F)
      if(refresh>0){cat(' . ')}
      
      if(smallSampleWarning){warning('Your strata sizes may be very small. Be careful interpreting stratified bootstrap outputs.')}
    }
    
    m_o_out <- m_o_out[2:nrow(m_o_out),] # in the case of bootstrapping, get rid of the overall fit in the "samples"
    
    attr(m_o_out,'indices') <- boot_sample_indices
    rm(boot_sample_indices)
  }
  
  m_brm <- match_brm_stan_optim(m_o_out=m_o_out,m_stan=m_stan,m_brm=m_brm,
                                m_o_overall = m_o_overall,
                                nSamples = nSamples,strata=strata
                                ,useBoot = useBoot)
  
  
  
  if(crossvals > 0){
    source('brm_model_oos_ll.R')
    
    smallSampleWarning <- F
    m_cv_all <- data.frame()
    
    for(. in 1:crossvals){ 
      
      m_cv_out <- brm_model_oos_ll(m_brm = m_brm
                                   , m_stan = m_stan
                                   # , m_o_overall = m_o_overall
                                   , strata = strata
                                   , nRuns = 2)
      
      m_cv_all <- rbind(m_cv_all,m_cv_out)
      
      rm(m_cv_out,m_brm_cv,loglik_oos,loglik_oos_dist,m_cv_2,m_cv_1)
      
      gc(verbose = F) ; Sys.sleep(.1) ; gc(verbose = F)
      
      if(refresh>0){cat(' . ')}
      
      if(smallSampleWarning){warning('Your strata sizes may be very small. Be careful interpreting cross-validation outputs.')}
    }
    
    m_brm$crossVals <- m_cv_all
    m_brm$criteria$logLik_oos <- m_cv_all$logLik_oos
    attr(m_brm$criteria$logLik_oos,'description') <- 
      'Out-of-sample log-likelihods. Mean and boostrapped -1SE and +1SE of pointwise [casewise] log likelihoods were calculated for each out-of-sample [test] set, then these values were multiplied by the size of the original data in order to extrapolate to the full-sample log-likelihoods.'
    m_brm$logLik_oos <- mean(m_cv_all$logLik_oos[,'50%'],trim = .25)
    attr(m_brm$logLik_oos,'description') <- paste('Robust central tendency measure of out-of-sample log likelihood. Of the',crossvals,'runs of random-split (90% train / 10% test) random cross-validation, the mean of the out-of-sample log-likelihoods was calculated after removing the highest and the lowest 25% of the distribution.')
    
  }else{
    m_brm$crossVals <- NA
    m_brm$logLik_oos <- NA
  }
  
  
  return(m_brm)
}
