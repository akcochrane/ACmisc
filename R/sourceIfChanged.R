

#' Source a file, but only if it changed from the last time you ran \emph{sourceIfChanged()} on the file
#'
#' Rather than sourcing your analysis script many times, which can be cumbersome, this script skips the
#' sourcing \emph{if} there is a \code{.log} file that matches the sourced file. If that matching \code{.log}
#' file exists, the corresponding \code{.RData} file is loaded. Otherwise, the file is sourced and both
#' \code{.log} and \code{.RData} files are created.
#'
#' @param filepath .R file to source
#'
#' The intended use is largely for when the user wants reproducible results (sometimes requiring
#' long analyses) while writing .Rmarkdown documents. While working with the main text, formatting,
#' and such, it's nice to just load pre-run data rather than waiting to source every time.
#'
#' There are important implications, such as the fact that only the \emph{sourced file} is checked for
#' changes, and therefore the old data will be loaded even if other things (e.g., global variables)
#' have been changed.
#'
#' @export
#'
sourceIfChanged <- function(filepath){

  ## TO DO: save the entire script as a variable as well

  filepathLength <- nchar(filepath)

  if( tolower(substr(filepath,filepathLength-1,filepathLength)) != tolower('.r')){
    stop('Please only try to source .R files')
  }

  file_contents <- readChar(filepath, file.info(filepath)$size)

  file_details <-   substr(filepath,1,filepathLength-2)

  logExists <- file.exists(paste0(file_details,'.log'))

  if(logExists){
    prev_log <- readChar(paste0(file_details,'.log'), file.info(paste0(file_details,'.log'))$size)

    if(prev_log != file_contents){
      logExists <- F
    }
  }

  if(!logExists){
    source(filepath)
    save.image(paste0(file_details,'.RData'))
    writeChar(file_contents,paste0(file_details,'.log'), nchar(file_contents))
  }
  load(paste0(file_details,'.RData'),envir = .GlobalEnv)
}
