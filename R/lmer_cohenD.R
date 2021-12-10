#' Find Cohen's D for lmer model
#' 
#' Calculate a heuristic Cohen's D for each fixed effect in an lmer model.
#' Uses a fairly naive calculation: \code{abs(b)/(SE*sqrt(df))}, where df is approximated
#' using a Kenward-Rogers approximation (package \code{pbkrtest} via \code{lmSupport::modelSummary}). 
#' Returns the \code{lmSupport::modelSummary} output with the addition of \code{$cohen_d}.
#' 
#' As with (& because of the underlying K-R approximation in) \code{lmSupport::modelSummary},
#' this can take a very long time for large or complex models.
#' 
#'
#' @param lmerMod Model fit by \code{lme4::lmer()}
#'
#' @export
#'
lmer_cohenD <- function(lmerMod){
  
  capture.output({suppressWarnings({
    m_summ <- lmSupport::modelSummary(lmerMod,Print = FALSE)
  })
                  }, file='NUL')
  
  m_cohenD <- abs(m_summ$KRAppox[,'Estimate'])/(m_summ$KRAppox[,'SE']*sqrt(m_summ$KRAppox[,'error df']))
  
  return(data.frame(m_summ$KRAppox,pseudo_cohenD = m_cohenD))
}