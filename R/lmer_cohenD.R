#' Find Cohen's D for lmer model
#' 
#' Calculate a heuristic Cohen's D for each fixed effect in an lmer model.
#' Uses a fairly naive calculation: \code{abs(b)/(SE*sqrt(df))}, where df is approximated
#' using a Kenward-Rogers approximation (package \code{pbkrtest} via \code{lmSupport::modelSummary}). 
#' Returns the \code{lmSupport::modelSummary} output with the addition of \code{$cohen_d}.
#'
#' @param lmerMod Model fit by \code{lme4::lmer()}
#'
#' @export
#'
lmer_cohenD <- function(lmerMod){
  
  m_summ <- lmSupport::modelSummary(lmerMod)
  
  m_cohenD <- abs(m_summ$KRAppox[,'Estimate'])/(m_summ$KRAppox[,'SE']*sqrt(m_summ$KRAppox[,'error df']))
  
  return(m_cohenD)
}