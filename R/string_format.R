##FUNCTIONS

#' Convert strings to lowercase in a smart way.
#'
#' @param x   A \code{character} string.
#' @param sep A \code{character} vector to specify multiple separator to remove
#'            during lowercase conversion (Default: sep = c(" ", "_")).
#' @return A \code{character} string converted in lowercase.
#' @author Yoann Pageaud.
#' @export
#' @examples smart.tolower(x = "A_QUIET_PLACE_PART_II")
#' @keywords internal

smart.tolower <- function(x, sep = c(" ", "_")){
  #Split words
  lapply(X = sep, FUN = function(i){
    x <<- unlist(strsplit(x, i))
  })
  #Check if a word is a roman number
  is.roman <- grepl(pattern = "^[IVMCXLD]{1}[IVMCXLD]+$", x = x)
  #Convert
  s <- vapply(
    X = seq_along(x), USE.NAMES = FALSE, FUN.VALUE = character(length = 1),
    FUN = function(i){
      if(!is.roman[i]){
        paste(substring(x[i], 1, 1), tolower(substring(x[i], 2)), sep = "")
      } else { x[i] }
    })
  paste(s, collapse = sep)
}
