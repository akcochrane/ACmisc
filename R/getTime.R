

#' Get the current date and time in an easy-to-use format
#' 
#' Returns a string with the date and time. 
#' 
#' Details
#' 
#' @param style style 1 uses letters to convey the month, style 2 uses numbers to convey the month.
#' @keywords  keywords
#' @export
#' @examples 
#' getTime()
#' 
getTime <- function(style=1){
  # style 1 uses month words, style 2 uses month numbers
  if (style==1){
    dateNow <- date()
    timeNow <- paste(substr(dateNow,21,24),substr(dateNow,5,7),substr(dateNow,9,10),'_',substr(dateNow,12,13),substr(dateNow,15,16),sep='')
  } else{
    dateNow <- Sys.time()
    timeNow <- paste(substr(dateNow,1,4),substr(dateNow,6,7),substr(dateNow,9,10),substr(dateNow,12,13),substr(dateNow,15,16),sep='')
  }
  return(timeNow)
  }