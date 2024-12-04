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
  if (is.numeric(q)) {
    q <- format(q, nsmall = n)
  }
  glue::glue("{s}={q}")
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

#' Construct middle part of filename of simulation run results
#'
#' @param alpha Real value less than 1, concentration parameter of the Pitman-Yor process
#' @param theta Real value, >= -alpha if alpha is between 0 (inclusive) and 1 (exclusive), or a positive integer multiple of |alpha| if alpha < 0
#' @param n Positive integer, number of nodes of the network i.e. order of corresponding graph
#' @param data.seed Positive integer, the seed for set.seed()
#' @param p Real value between 0 and 1 inclusive, prior probability of finite regime
#' @export
filename_body <- function(alpha, theta, n, data.seed, p) {
  str0 <- append_format("alpha", alpha, 2)
  str1 <- append_format("theta", theta, 2)
  str2 <- append_format("n", n, 0)
  str3 <- append_format("data.seed", data.seed, 0)
  str4 <- append_format("p", p, 2)
  paste("sim_dc", str0, str1, str2, str3, str4, sep = "_")
}
