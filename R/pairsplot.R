
#' Create a \code{psych::pairs.panels()} plot with standard formatting
#'
#' The default title includes approximate Bayes Factors-based estimates of regions providing "At least moderate
#' evidence for the lack of a correlation" or "at least moderate evidence for a correlation." The BF approximation
#' uses the \code{BayesFactor::ttest.tstat()} function, and finds the r value associated with a T value associated with a
#' raw BF <.33333 or >3. \strong{Beware}: The \emph{n} used in these calculations is the number of rows of the data frame.
#' If you have missing data, these estimates could be wildly off.
#'
#' @param d Data frame to be plotted
#' @param method Correlation method to be used: product-moment ('pearson') or rank ('spearman') ?
#' @param main Plot title, as in plot()
#' @param stars Do you want magic stars? (same as \code{pairs.panels()})
#' @param color What color do you want in the histogram bars?
#'
#' @export
#'
#' @examples
#' d <- iris[1:100,c('Sepal.Length','Sepal.Width','Petal.Length','Petal.Width')]
#' d$isoDist <- ACmisc::d1iso(d)
#' pairsplot(d)
#' 
pairsplot <- function(d, method = 'Spearman', main = '',stars=T,color='darkgrey'){

  library(psych)

  par(family='serif')

  if(nchar(main)==0){

    alt_thresh <- suppressMessages({suppressWarnings({
      optimise(
      function(r_val){
        (3-BayesFactor::ttest.tstat(psych::r2t(r_val,nrow(d)),nrow(d),simple = T) )^2
      },c(.000001,.99999))$minimum
  })})

    null_thresh <- suppressMessages({suppressWarnings({
    optimise(
      function(r_val){
        (.333333333333-BayesFactor::ttest.tstat(psych::r2t(r_val,nrow(d)),nrow(d),simple = T) )^2
      },c(.000001,.99999))$minimum
    })})

    main = paste0('Pairs Plot, ',method,'\nApproximate null acceptance region |r|<',round(null_thresh,2),'; null rejection region |r|>',round(alt_thresh,2))
    }

  for(curVar in colnames(d)){
    if(is.character(d[,curVar])){d[,curVar] <- NULL}
    if(is.factor(d[,curVar])){
      if(length(unique(d[,curVar])) > 2){d[,curVar] <- NULL}
    }
  }

  pairs.panels(d,
               method = tolower(method),
               ci=T,
               lm = T,
               stars = stars,
               ellipses = F,
               alpha = .045,
               hist.col = color,
               # breaks = round(nrow(na.omit(d))*.3333),
               main = main)
}
