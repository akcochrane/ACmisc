

#' Extract first dimension using nonmetric scaling
#'
#' Implements \code{\link[MASS]{isoMDS}(dist(data,dist_method),k=1)}, with
#' attributes of correlations with (scaled) rowMeans and with PCA component-1.
#' 
#' Will break if identical cases are present.
#'
#' @param data Must be numeric
#' @param dist_method 'manhattan' (default; L1 norm) or 'euclidean' (AKA OLS or L2 norm).
#'
#' @export
#' 
#' @examples 
#' d <- iris[1:100,c('Sepal.Length','Sepal.Width','Petal.Length','Petal.Width')]
#' d$isoDist <- d1iso(d)
#' ACmisc::pairsplot(d)
#'
d1iso <- function(data,dist_method = c('manhattan','euclidean')){

  data_dist <- dist(data,dist_method[1])

  dat_d1 <- MASS::isoMDS(data_dist,k=1)$points

  dat_d1 <- (dat_d1 - median(dat_d1,na.rm=T))/mad(dat_d1,na.rm=T)

  scaleRowMeans <- rowMeans(scale(data))

  if(cor(dat_d1,scaleRowMeans,use='complete') < 0){dat_d1 <- -dat_d1}

  attr(dat_d1,'pearson_w_rowMean') <- cor(dat_d1,scaleRowMeans,use='complete')
  attr(dat_d1,'pearson_w_component') <- cor(dat_d1,psych::principal(data,nfactors = 1)$scores,use='complete')

 return(dat_d1)
}
