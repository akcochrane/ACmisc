
#' Use delta-BIC to approximate Bayes Factors for model comparison
#'
#' Calculates BIC for the entire model, and for the entire set of models
#' dropping each fixed effect one at a time. Uses the change in BIC between
#' models to approximate Bayes Factors (essentially penalized \emph{Likelihood Ratios}, or
#' evidence for one model over another). This approach is promoted by, for example,
#' Kass and Raftery (1995; 10.1080/01621459.1995.10476572) and Wagenmakers (2007; 10.3758/BF03194105).
#'
#' NOT intended for relations between two variables (in that case, use \code{robustCor} or similar).
#'
#' Could use further testing to ensure that formula extraction is correct.
#'
#' Follows the recommendations of \emph{Wagenmakers (2007; PBR)}, e.g., for LMEMs calculates BIC
#' using \emph{-2*logLik(modIn) + attr(logLik(modIn),'df')*log(sum(ngrps(modIn)))}
#'
#' @note
#' BIC_dropped of the Intercept term is the \strong{full model BIC}. BIC_dropped of all
#'  other coefficients are the BIC of the model without that coefficient.
#'  See comment(output) for more information.
#'  
#'  \strong{Known Bug} preventing use with non-multi-level models. In other words,
#'  it is only currently functional for multi-level models (i.e., \code{lmer()} or \code{glmer()})
#'
#' @param modIn \code{lm}, \code{glm}, \code{rlm}, \code{lmer} or \code{glmer} model.
#' @param logBase Base of the logBF
#'
#' @export
#'
#' @examples
#' d <- data.frame(y = rnorm(60),grup = sample(c('a','b','c','d'),60,replace = T),arb = sample(c('e','f','g','h'),60,replace = T))
#' d$x  <- d$y + rnorm(60,0,.2)
#'
#' m_lmer <- lme4::lmer(y~x*arb+(1|grup),d,REML = F)
#' BICBF(m_lmer)
#'
#' m_glmer <- lme4::glmer(y~x*arb+(1|grup),d,family=gaussian())
#' BICBF(m_glmer)
#'
#' m_lm <- lm(ravens ~ alert * corsi + isAdult, dat_cochraneEtAl_2019_PLOSOne)
#' BICBF(m_lm)
#'
#' m_rlm <- MASS::rlm(ravens ~ alert * corsi + isAdult, dat_cochraneEtAl_2019_PLOSOne)
#' BICBF(m_rlm)
#'
#' m_glm <- glm(y ~ x*arb,d,family=gaussian())
#' BICBF(m_glm)
#' 
BICBF <- function(modIn,logBase=3){

  require(MASS)

  logVarName <- paste0('BFlog',round(logBase,2))

  if(class(modIn)[1] == 'lmerMod' || class(modIn)[1] == 'glmerMod'){ ## lme4 mixed-effects models (i.e., (g)lmer models)

    require(lme4)
    ## convert the concise formula into a [factor-wise] verbose one:

    dat <- as.data.frame(modIn@pp$X)
    colnames(dat) <- colnames(modIn@pp$X)
    dat[,'(Intercept)'] <- NULL


    termLabels <- attr(terms(formula(modIn)),'term.labels')
    ranEfs <- paste('(',termLabels[grep('|',termLabels,fixed=T)],')')

    modForm <- paste(formula(modIn)[2],'~',
                     paste(c(colnames(dat),ranEfs),
                           collapse=' + '))

    dat <- data.frame(dat,modIn@frame[,setdiff(colnames(modIn@frame),colnames(dat))])


    if(family(modIn)$family=='gaussian'){
      modInternal <- lmer(modForm,dat,REML = F)
    }else{
      modInternal <- glmer(modForm,dat,family(modIn))
    }

    modIn <- modInternal

    ## Move on to the actual fitting / comparisons
    outTab <- data.frame(summary(modIn)$coefficients)
    outTab$BIC_dropped <- NA
    outTab[,logVarName] <- NA

    for(curPred in rownames(outTab)){
      if(curPred == '(Intercept)'){
        outTab[curPred,'BIC_dropped'] <- -2*logLik(modIn) + attr(logLik(modIn),'df')*log(sum(ngrps(modIn)))
      }else{

        suppressMessages({
          suppressWarnings({
            curMod <- update(modIn, as.formula(paste(' . ~ . -', curPred) ))
            outTab[curPred,'BIC_dropped'] <- -2*logLik(curMod)+
              + attr(logLik(curMod),'df')*log(sum(ngrps(curMod)))
          })})

        # https://statswithr.github.io/book/bayesian-model-selection.html

        outTab[curPred,logVarName] <- log(
          exp((outTab[curPred,'BIC_dropped'] - outTab['(Intercept)','BIC_dropped'])/2)
          ,base=logBase
        )

        # rm(modDiff)
      }
    }

  }else{ # before here is the (g)lmer option

    ### this is a total hack, but extracting this model frame seems easier from a rlm model than from a (g)lm model
    if(class(modIn)[1] == 'lm' || class(modIn)[1] == 'glm'){
      dat <- data.frame(modIn$model)
      suppressWarnings({suppressMessages({dat <- data.frame(
        rlm(formula(modIn),data=dat,maxit=5)$x
      )})})

    }
    if(class(modIn)[1] == 'rlm'){ dat <- data.frame(modIn$x) }

    colnames(dat) <- names(coef(modIn))

    ## THIS SECTION MAKES IT BREAK WHEN THERE'S ONLY ONE PREDICTOR

    dat[,grep('Intercept',colnames(dat))] <- NULL # remove intercept offset from model frame

    # termLabels <- attr(terms(formula(modIn)),'term.labels') ## unnecessary?

    modForm <- paste(formula(modIn)[2],'~',
                     paste(colnames(dat),
                           collapse=' + '))

    dat <- data.frame(dat,modIn$model[,setdiff(colnames(modIn$model),colnames(dat))])
    colnames(dat)[ncol(dat)] <- setdiff(colnames(modIn$model),colnames(dat))

    ## this is sloppy. should use switch
    if(class(modIn)[1] == 'rlm'){
      suppressMessages({ suppressWarnings({
        modInternal <- rlm(as.formula(modForm),dat,k2=modIn$k2,weights = modIn$weights,maxit=200)
      })})}
    if(class(modIn)[1] == 'glm'){
      modInternal <- glm(as.formula(modForm),data=dat,family = modIn$family,weights = modIn$weights)
    }
    if(class(modIn)[1] == 'lm'){
      modInternal <- lm(modForm,data=dat)
    }


    ## Move on to the actual fitting / comparisons
    outTab <- data.frame(summary(modInternal)$coefficients[,1:3])
    outTab$BIC_dropped <- NA
    outTab[,logVarName] <- NA
    #

    for(curPred in rownames(outTab)){
      if(curPred == '(Intercept)'){
        outTab[curPred,'BIC_dropped'] <- BIC(modInternal)
      }else{

        # suppressMessages({
        # suppressWarnings({
        curMod <- update(modInternal, as.formula(paste(' . ~ . -', curPred) ))
        outTab[curPred,'BIC_dropped'] <- BIC(curMod)
        # })})

        # https://statswithr.github.io/book/bayesian-model-selection.html

        outTab[curPred,logVarName] <- log(
          exp((outTab[curPred,'BIC_dropped'] - outTab['(Intercept)','BIC_dropped'])/2)
          ,base=logBase
        )

      }
    }

  } ## close the rlm/glm/lm options

  comment(outTab) <- '"BIC_dropped" associated with Intercept is full model BIC.
  Bayes factor calculation is approximated as exp((BIC1-BIC2)/2).'

  outTab[,c(logVarName,'BIC_dropped')] <- round(outTab[,c(logVarName,'BIC_dropped')],3)

  outTab[,grep('value',colnames(outTab))] <- round(outTab[,grep('value',colnames(outTab))],3)

  return(outTab)

}
