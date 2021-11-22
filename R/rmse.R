#' Find Root Mean Squared Error
#' 
#' Simply finds the RMSE between vect1 and vect2 (of the same length). 
#' Often, when assessing fit values, this is interpreted as the standard 
#' deviation of fits (e.g., vect2) around true values (e.g., vect1).
#'
#' @param vect1 A numeric vector
#' @param vect2 Another numeric vector, of the same length as \code{vect1}
#'
#' @export
#'
rmse <- function(vect1,vect2){
  
  return(
    sqrt(mean((vect1-vect2)^2,na.rm=T))
  )
  
}
