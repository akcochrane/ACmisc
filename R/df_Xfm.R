#' Transform a numeric data frame or matrix using the Yeo-Johnson algorithm
#'
#' Returns a data frame, along with many supplementary attributes
#' (e.g., correlations, transformation parameters, simulate from correlation matrix).
#'
#' @param x Numeric data frame or matrix
#' @param keepMedian	Should the result be shifted so that it has the same median as the input vector?
#' @param keepMAD Should the result be scaled so that it has the same dispersion (median absolute deviation) as the input vector?
#'
#' @export
#'
#' @examples
#' x <- data.frame(w = 'char',x=rnorm(20),y=rgamma(20,3,2),z=c(rnorm(19,3,3),NA))
#' x_transformed <- df_xfm(x)
#' head(x_transformed)
#' attributes(x_transformed)
#'
df_xfm <- function(x, keepMedian = T, keepMAD = T){

  lambda <- c()
  df_Xfmd <- data.frame(matrix(NA,nrow(x),0))
  rownames(df_Xfmd) <- rownames(x)

  num_cols <- rep(F,ncol(x))

  for(curCol in 1:ncol(x)){
    try({
      x_Xfmd <- ACmisc::YeoJohn(x[,curCol], keepMedian=keepMedian, keepMAD=keepMAD)
      lambda <- c(lambda,attr(x_Xfmd,'lambda'))
      df_Xfmd <- cbind(df_Xfmd,x_Xfmd)
      num_cols[curCol] <- T
    },silent = T)

    if(ncol(df_Xfmd) < curCol){df_Xfmd <- cbind(df_Xfmd,x[,curCol])}

    colnames(df_Xfmd)[ncol(df_Xfmd)] <- colnames(x)[curCol]
  }

  attr(df_Xfmd,'var_yeoJohnson_lambda') <- lambda
  attr(df_Xfmd,'var_pearson') <- cor(df_Xfmd[,num_cols],method = 'pearson',use='pairwise')
  attr(df_Xfmd,'var_spearman') <- cor(df_Xfmd[,num_cols],method = 'spearman',use='pairwise')
  attr(df_Xfmd,'var_median') <- apply(df_Xfmd[,num_cols],2,median,na.rm=T)
  attr(df_Xfmd,'var_MAD') <- apply(df_Xfmd[,num_cols],2,mad,na.rm=T)

  attr(df_Xfmd[,num_cols],'orig_mean') <- apply(x[,num_cols],2,mean,na.rm=T)
  attr(df_Xfmd[,num_cols],'orig_SD') <- apply(x[,num_cols],2,sd,na.rm=T)
  attr(df_Xfmd[,num_cols],'orig_skew') <- apply(x[,num_cols],2,skew,na.rm=T)
  attr(df_Xfmd[,num_cols],'orig_kurtosis') <- apply(x[,num_cols],2,kurtosi,na.rm=T)

  attr(df_Xfmd,'numeric_cols') <- num_cols

  attr(df_Xfmd,'sim_function') <- function(nSims,df_Xfmd){
    simDat <- MASS::mvrnorm(nSims,
                            eval(attr(df_Xfmd,'var_median')),
                            cor2cov(attr(df_Xfmd,'var_pearson'),
                                    attr(df_Xfmd,'var_MAD')))
    colnames(simDat) <- colnames(df_Xfmd)[attr(df_Xfmd,'numeric_cols')]
    return(data.frame(simDat))
  }

  return(df_Xfmd)

}
