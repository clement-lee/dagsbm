#' Post-MCMC process results of fitting DAG-SBM
#'
#' @param obj List, object returned from \code{\link{mcmc_dagsbm}}
#' @param proc_time A proc_time object, usually from \code{\link{system.time}}
#' @importFrom tibble as_tibble
#' @importFrom stringr str_c
#' @importFrom utils sessionInfo
#' @export
post_mcmc_process <- function(obj, proc_time) {
  iter <- obj$scalars$iter
  thin <- obj$scalars$thin
  burn <- obj$scalars$burn
  obj$scalars$seconds <- proc_time["elapsed"]
  obj$scalars$spi <- proc_time["elapsed"] / (iter * thin + burn)
  names_tbl <- function(x, prefix, make_integer) {
    ## assign colnames & make tibble, optionally convert to integer
    n <- ncol(x)
    if (make_integer) {
      x <- as.data.frame(apply(x, 2, as.integer))
    }
    colnames(x) <- stringr::str_c(prefix, seq(n))
    x <- tibble::as_tibble(x)
    x
  }
  obj$scalars <- obj$scalars |> tibble::as_tibble()
  obj$pars <- obj$pars |> tibble::as_tibble()
  obj$point_est <- obj$point_est |> tibble::as_tibble()
  obj$Z <- obj$Z |> names_tbl("Z_", make_integer = TRUE)
  obj$phi <- obj$phi |> names_tbl("phi_", make_integer = TRUE)
  obj$xi <- obj$xi |> names_tbl("xi_", make_integer = FALSE)
  ## obj$A & obj$Y stays unconverted
  obj$sessionInfo <- utils::sessionInfo()
  obj
}

#' Read results of fitting an SBM to a set of simulated data
#'
#' @param alpha Real value less than 1, concentration parameter of the Pitman-Yor process
#' @param theta Real value, >= -alpha if alpha is between 0 (inclusive) and 1 (exclusive), or a positive integer multiple of |alpha| if alpha < 0
#' @param n Positive integer, number of nodes of the network i.e. order of corresponding graph
#' @param p Real value between 0 and 1 inclusive, prior probability of finite regime
#' @param data.seed Positive integer, the seed for set.seed()
#' @param dir Boolean, is the SBM for directed graphs (TRUE) or directed acyclic graphs (FALSE)?
#' @importFrom glue glue
#' @importFrom here here
#' @importFrom readr read_rds
#' @seealso \code{\link{obtain_spi}}, \code{\link{compute_bf}} and \code{\link{plot_topo_wrapper}} with the same set of arguments as they call this function
#' @export
read_results <- function(alpha, theta, n, p, data.seed, dir) {
  str0 <- filename_body(alpha, theta, n, data.seed, p) # diff arg ordering
  str0 <- glue::glue(".*{str0}_.*{ifelse(dir,'dir','dag')}\\.rds")
  v0 <- list.files("results/", str0)
  if (length(v0) == 0L) {
    stop("read_results: there's no file with the provided combination of parameter values.")
  } else if (length(v0) > 1L) {
    print(v0)
    message("read_results: there's more than 1 file with the provided combination of parameter values. The chronologically earlier one will be used.")
    v0 <- v0[1L]
  }
  here::here("results", v0) |> readr::read_rds()
}

#' Restructure data frame of MCMC samples for pivot_topo() and subsequent plotting
#'
#' @param df Data frame of MCMC samples, with at least the following column names: K, k, gamma, theta, alpha, a0, b0, finite
#' @param var_names Character vector with one or more of the following values (in quotes): Kn, k, gamma, theta, alpha, log(gamma), a, b
#' @importFrom dplyr transmute filter mutate
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#' @importFrom forcats as_factor
#' @export
pivot_pars <- function(df, var_names) {
  ## var_names is after renaming {a0,b0} to {a,b}
  df |>
  dplyr::transmute( # order SHOULD be preserved
    Kn = .data$K,
    k = ifelse(.data$finite == 0.0, as.integer(NA), .data$k),
    gamma = ifelse(.data$finite == 0.0, as.numeric(NA), .data$gamma),
    theta = ifelse(.data$finite == 1.0, as.numeric(NA), .data$theta),
    alpha = ifelse(.data$finite == 1.0, as.numeric(NA), .data$alpha),
    `log(gamma)` = log(.data$gamma),
    a = .data$a0,
    b = .data$b0
  ) |>
  tidyr::pivot_longer(.data$Kn:.data$b) |>
  dplyr::filter(.data$name %in% var_names) |>
  dplyr::mutate(name = forcats::as_factor(.data$name))
}

#' Restructure data frame for topological ordering plotting
#'
#' @param df Data frame, normally output of pivot_pars
#' @param v Integer vector, a permutation of 1, 2, ..., length(v)
#' @importFrom tibble tibble
#' @importFrom glue glue
#' @importFrom dplyr sample_n everything group_by summarise count mutate ungroup left_join
#' @importFrom tidyr pivot_longer
#' @export
pivot_topo <- function(df, v = NULL) {
  ## df ~ that for pivot_pars()
  if (!(missing(v) || is.null(v))) {
    if (!all.equal(sort(v), seq_along(v))) {
      stop("pivot_topo: when v is non-null, it has to be a permutation of {1,...,n}, where n is the number of nodes.")
    } else {
      df1 <- tibble::tibble(name = glue::glue("phi_{seq_along(v)}"), position = as.integer(v))
    }
  } else {
    v <- df |> dplyr::sample_n(1) |> unlist()
    df1 <- tibble::tibble(name = names(v), position = as.integer(v))
  }
  df2 <-
    df |>
    tidyr::pivot_longer(dplyr::everything()) |>
    dplyr::group_by(.data$name) |>
    dplyr::summarise(mean = mean(.data$value))
  df |>
    tidyr::pivot_longer(dplyr::everything()) |>
    dplyr::count(.data$name, .data$value) |>
    dplyr::group_by(.data$name) |>
    dplyr::mutate(density = .data$n / sum(.data$n)) |>
    dplyr::ungroup() |>
    dplyr::left_join(df1, "name") |>
    dplyr::left_join(df2, "name")
}

#' Calculate Bayes factor of finite regime to infinite regime of the DAG-SBM, for a set of simulated data
#'
#' @param alpha Real value less than 1, concentration parameter of the Pitman-Yor process
#' @param theta Real value, >= -alpha if alpha is between 0 (inclusive) and 1 (exclusive), or a positive integer multiple of |alpha| if alpha < 0
#' @param n Positive integer, number of nodes of the network i.e. order of corresponding graph
#' @param p Real value between 0 and 1 inclusive, prior probability of finite regime
#' @param data.seed Positive integer, the seed for set.seed()
#' @param dir Boolean, is the SBM for directed graphs (TRUE) or directed acyclic graphs (FALSE, default)?
#' @importFrom stringr str_c
#' @seealso \code{\link{obtain_spi}}, \code{\link{read_results}} and \code{\link{plot_topo_wrapper}} with the same set of arguments
#' @export
compute_bf <- function(alpha, theta, n, p, data.seed, dir = FALSE) {
  obj0 <- read_results(alpha, theta, n, p, data.seed, dir)
  bf <- odds(obj0$pars$finite |> mean()) / odds(obj0$scalars$p_fin)
  A <- floor(log10(bf))
  if (A >= -3.0) {
    bf <- format(bf, digits = 2, nsmall = 2)
  } else {
    bf <- stringr::str_c(format(bf / 10 ^ A, digits = 2, nsmall = 2), "\\times10^{", A, "}")
  }
  stringr::str_c("$", bf, "$")
}

#' Obtain seconds per iteration of MCMC runs
#'
#' @param alpha Real value less than 1, concentration parameter of the Pitman-Yor process
#' @param theta Real value, >= -alpha if alpha is between 0 (inclusive) and 1 (exclusive), or a positive integer multiple of |alpha| if alpha < 0
#' @param n Positive integer, number of nodes of the network i.e. order of corresponding graph
#' @param p Real value between 0 and 1 inclusive, prior probability of finite regime
#' @param data.seed Positive integer, the seed for set.seed()
#' @param dir Boolean, is the SBM for directed graphs (TRUE) or directed acyclic graphs (FALSE)?
#' @seealso \code{\link{compute_bf}}, \code{\link{read_results}} and \code{\link{plot_topo_wrapper}} with the same set of arguments
#' @export
obtain_spi <- function(alpha, theta, n, p, data.seed, dir) {
  signif(read_results(alpha, theta, n, p, data.seed, dir)$scalars$spi, 2)
}
