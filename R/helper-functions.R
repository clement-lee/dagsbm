#' Calculate odds from probability
#' 
#' @param x Real value between 0 and 1 exclusive
#' @export
odds <- function(x) x / (1.0 - x)

#' Format string for filename purposes
#'
#' @param s String, name of quantity
#' @param q Numeric, actual value of quantity
#' @param n Positive integer, number of decimal places desired
#' @importFrom glue glue
#' @export
append_format <- function(s, q, n) {
  glue::glue("{s}={format(q, nsmall = n)}")
}

#' Check if object exists, with error message if not
#'
#' @param string Scalar string
#' @importFrom glue glue
#' @export
stop_wrapper <- function(string) {
  if (!exists(string)) {
    stop(glue::glue("{string} has to be supplied"))
  }
}

#' Remove dash in current date & make string
#'
#' @importFrom stringr str_remove_all
#' @export
yyyymmdd <- function() {
  stringr::str_remove_all(as.character(Sys.Date()), "-")
}
