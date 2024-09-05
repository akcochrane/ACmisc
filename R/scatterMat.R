#' Create a matrix of scatter plots, with each row and each variable containing one variable
#'
#' @param df Data frame
#' @param colNames Character vector of variable names to create the columns of the correlation matrix
#' @param rowNames Character vector of variable names to create the rows of the correlation matrix
#' @param corMethod Should the reported correlation be "pearson" or "spearman"?
#'
#' @export
#'
#' @examples
#' 
#' scatterMat(iris, colNames = c('Sepal.Width','Petal.Width'),c('Sepal.Length','Petal.Length'))
#' 
scatterMat <- function(df, colNames, rowNames, corMethod = 'spearman'){
  # df <- d_subj[d_subj$Gender < 3 , ]
  library(ggplot2)
  library(gridExtra)
  
  gList <- list();  for(curRow in rowNames){
    for(curCol in colNames){
      
      gList[[paste0(curRow,'_',curCol)]] <-
        ggplot(df, aes_string(y = curCol, x = curRow)) +
        theme_bw() +
        geom_smooth(method='lm') +
        geom_point() +
        labs(x = paste0(curRow , '\nr =' ,round(cor(df[,curRow],df[,curCol], method = corMethod, use = 'complete'),3))) +
        theme(
          text = element_text(family = 'serif')
          ,axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10) )
    }
    
  }
  grid.arrange(grobs=gList , ncol = length(colNames))
}