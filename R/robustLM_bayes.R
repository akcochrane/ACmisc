#' Fit a robust model, with fully Bayesian model comparisons
#'
#' Fits a median linear regression, with random effects (i.e., mixed-effects model) if desired,
#' using \code{\link[brms]{brm}} and the asymmetric laplace distribution (fitting quantile = .5, i.e., median).
#' Then fits a set of models dropping each fixed-effect term in turn, and uses
#' bridge sampling to estimate the Bayes Factors of the full model relative to each
#' sub-model. Returns the full model, with a data frame \code{model$bayes_factors} with
#' each fixed effect's log (base 3) Bayes Factor.
#'
#' Bridge sampling is run 5 times for each effect, and the most conservative (i.e., equivocal)
#' of the resulting BF is returned.
#'
#' While using \code{algorighm = 'fullrank'} speeds up the model (over the default \code{algorithm  = 'sampling'} by a factor of over 3.5 in preliminary tests), it would be advised
#' to fit and compare several fitting runs to assess the robustness of the resulting coefficients and Bayes Factors. In brief tests, some
#' values have been quite consistent across fitting runs and others have not when using \code{fullrank} or \code{meanfield}; in general, waiting for
#' \code{sampling} provides the most robust estimates.
#'
#' @seealso
#' \code{\link[brms]{brm}}, \code{\link[brms]{bayes_factor}}, and \code{\link[brms]{brmsfamily}}.
#'
#' @param formula Formula, as in \code{lm()}, \code{brm()} or \code{lmer()}
#' @param data Data, as in \code{lm()} or \code{lmer()}
#' @param ...  Additional arguments to pass to \code{\link[brms]{brm}}. For example, if multiple cores are available, \code{cores=2} or more should drastically speed up performance.
#' @param iter Number of iterations. More is better, within reasonable time constraints
#' @param showProgress Show a small progress bar?
#' @param quiet Logical. To show model progress, or hide it.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This fairly simple example takes about 160 seconds on a good computer.
#' m <- robustLM_bayes(mot ~ ufov * enum, dat_cochraneEtAl_2019_PLOSOne) 
#'
#' m$bayes_factors
#' }
robustLM_bayes <- function(formula,data,...,iter = 10000 , showProgress = T, quiet = T){

  ## test with MEM; figure out how to not have grouping variables in the BF comparisons

  formIn <- formula ; rm(formula)
  dataIn <- data ; rm(data)
  require(brms)
  if(showProgress){pb <- txtProgressBar(max=3,width = 3,style=3)}
  suppressWarnings({ suppressMessages({
    m_full <- brm(brmsformula(formIn,quantile=.5,family=asym_laplace())
                  ,dataIn
                  ,iter = iter
                  ,silent = quiet
                  ,refresh = as.numeric(!quiet)
                  ,save_pars = save_pars(all=T)
                  # ,save_all_pars = T
                  , ...)
  })})
  if(showProgress){setTxtProgressBar(pb, 1)}
  m_dropOne <- list()
  bayes_factors <- data.frame()

  fixefNames <- attr(attr(m_full$data,'terms'),'term.labels')
  fixefNames <- fixefNames[2:length(fixefNames)]

  halfFixEf <- fixefNames[round(length(fixefNames) / 2)]

  for(curDrop in fixefNames){


    suppressWarnings({ suppressMessages({
      m_dropOne[[curDrop]] <- update(m_full
                                     , formula = as.formula(paste(' . ~ . -', curDrop) )
                                     ,silent = quiet
                                     ,refresh = as.numeric(!quiet)
                                     ,save_pars = save_pars(all=T)
                                     # ,save_all_pars = T
      )
    }) })


    bf_tmp <- c() ; for(curBF in 1:5){ # this re-runs the bridge sampling several times, and returns the most conservative (equivocal) Bayes Factor
      bf_tmp[curBF] <- NA
      try({
      bf_tmp[curBF] <- log(bayes_factor(m_full
                                        ,m_dropOne[[curDrop]]
                                        ,silent = quiet
      )$bf,base=3)
      })
    }
    if(all(is.na(bf_tmp))){stop('No bridge sampling runs succeeded.')}
    if(showProgress){ if(curDrop == halfFixEf){setTxtProgressBar(pb, 2)} }
    # print(bf_tmp) ## for testing
    bf_tmp <- na.omit(bf_tmp)

    bayes_factors <- rbind(bayes_factors
                           ,data.frame(predictor=curDrop
                                       ,BFlog3 = bf_tmp[which.min(abs(bf_tmp))])
    )
  }

  m_full$dropOneMods <- m_dropOne

  if(showProgress){ setTxtProgressBar(pb, 3) }

  
  if(showProgress){ close(pb) }
  m_full$bayes_factors <- bayes_factors
  return(m_full)
}
