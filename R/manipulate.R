#' Return a sorted vector of nodes id
#'
#' @param g An igraph object of a DAG
#' @param random Boolean, whether the order of selected nodes is randomised in the process
#' @return A data frame with two columns: "id" is the names of nodes in g, and "id_num" is the topological ordering
#' @examples
#' df0 <- data.frame(from = c("a", "b"), to = c("b", "c"), stringsAsFactors = FALSE)
#' g0 <- igraph::graph_from_data_frame(df0, directed = TRUE)
#' topo_sort_kahn(g0)
#' @importFrom igraph as_data_frame
#' @importFrom igraph V
#' @importFrom igraph as_adjacency_matrix
#' @importFrom igraph is_dag
#' @export
topo_sort_kahn <- function(g, random = FALSE) {
  if (!is_dag(g)) {
    stop("g has to be a DAG")
  }
  e0 <- igraph::as_data_frame(g)
  names(e0) <- c("citing", "cited")
  v <- names(igraph::V(g))
  v1 <- sort(unique(c(e0$citing, e0$cited)))
  l <- setdiff(v, v1)
  s0 <- sort(setdiff(e0$citing, e0$cited))
  while (length(s0) > 0L) {
    if (random) {
      n0 <- sample(s0, 1L)
      s0 <- s0[s0 != n0]
    } else {
      n0 <- s0[1L]
      s0 <- s0[-1L]
    }
    l <- c(l, n0)
    ## outgoing edges of n0
    i0 <- e0$citing == n0
    e1 <- e0[i0, , drop = FALSE]
    e0 <- e0[!i0, , drop = FALSE]
    if (nrow(e1) != 0L) {
      e2 <- setdiff(e1$cited, e0$cited)
      if (random) {
        e2 <- sample(e2, length(e2))
      }
      s0 <- c(s0, e2)
    }
  }
  if (nrow(e0) > 0L) {
    stop("topo_sort_kahn: graph has at least 1 cycle")
  }
  o <- match(l, v)
  a <- as_adjacency_matrix(g)[o, o]
  data.frame(
    id = l,
    id_num = seq_along(l),
    stringsAsFactors = FALSE
  )
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
#' @seealso \code{\link{obtain_spi}}, \code{\link{compute_bf}} and \code{\link{plot_topo_wrapper}} with the same set of arguments
#' @export
read_results <- function(alpha, theta, n, p, data.seed, dir) {
  alpha0 <- format(alpha, nsmall = 2)
  theta0 <- format(theta, nsmall = 2)
  if (is.numeric(p)) {
    p <- format(p, nsmall = 2)
  }
  string0 <- glue::glue(".*sim_dc_alpha={alpha0}_theta={theta0}_n={n}_data.seed={data.seed}_p={p}_.*{ifelse(dir,'dir','dag')}\\.rds")
  v0 <- list.files("results/", string0)
  if (length(v0) == 0L) {
    stop("read_results: there's no file with the provided combination of parameter values.")
  } else if (length(v0) > 1L) {
    print(v0)
    stop("read_results: there's more than 1 file with the provided combination of parameter values. Check if there's anything wrong.")
  }
  here::here("results", v0) |> readr::read_rds()
}

#' Obtain point estimate of the allocation vector
#'
#' @param Z.mat Data frame of MCMC samples of the node memberships
#' @param g igraph object of the underlying graph
#' @param ... other arguments passed to salso::salso()
#' @importFrom salso salso
#' @importFrom mcclust comp.psm
#' @importFrom igraph as_adjacency_matrix
#' @export
obtain_point_est <- function(Z.mat, g, ...) {
  ## Z.mat is $Z of an MCMC output
  Z.mat <- as.matrix(Z.mat)
  print(Sys.time())
  print("calculating point estimate")
  est <- salso::salso(Z.mat, ...)
  com <- list(membership = est) # for network plot
  class(com) <- "communities"
  cN <-
    est |>
    table() |>
    cumsum() # for plot_adj_psm()
  ord <-
    est |>
    rank(ties.method = "random") |>
    order()
  print(Sys.time())
  print("calculating co-clustering matrix")
  psm <-
    Z.mat |>
    mcclust::comp.psm() |>
    reorder_dense(ord-1L) |> # 0-indexing
    as.data.frame.table()
  print(Sys.time())
  print("calculating adjacency matrix")
  adj <-
    g |>
    igraph::as_adjacency_matrix() |>
    reorder_sparse(ord-1L) |> # 0-indexing
    as.matrix() |>
    as.data.frame.table()
  list(
    est = est,
    com = com,
    cN = cN,
    psm = psm,
    adj = adj
  )
}

#' Calculate odds from probability
#' 
#' @param x Real value between 0 and 1 exclusive
#' @export
odds <- function(x) x / (1.0 - x)

#' Calculate Bayes factor of finite regime to infinite regime of the DAG-SBM, for a set of simulated data
#'
#' @param alpha Real value less than 1, concentration parameter of the Pitman-Yor process
#' @param theta Real value, >= -alpha if alpha is between 0 (inclusive) and 1 (exclusive), or a positive integer multiple of |alpha| if alpha < 0
#' @param n Positive integer, number of nodes of the network i.e. order of corresponding graph
#' @param p Real value between 0 and 1 inclusive, prior probability of finite regime
#' @param data.seed Positive integer, the seed for set.seed()
#' @param dir Boolean, is the SBM for directed graphs (TRUE) or directed acyclic graphs (FALSE, default)?
#' @seealso \code{\link{obtain_spi}}, \code{\link{read_results}} and \code{\link{plot_topo_wrapper}} with the same set of arguments
#' @export
compute_bf <- function(alpha, theta, n, p, data.seed, dir = FALSE) {
  obj0 <- read_results(alpha, theta, n, p, data.seed, dir)
  bf <- odds(obj0$pars$finite |> mean()) / odds(obj0$scalars$p_fin)
  A <- floor(log10(bf))
  if (A >= -3.0) {
    bf <- paste0("$", format(bf, digits = 2, nsmall = 2), "$")
  } else {
    bf <- paste0("$", format(bf / 10 ^ A, digits = 2, nsmall = 2), "\\times10^{", A, "}$")
  }
  bf
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
obtain_spi <- function(alpha, theta, n, p, data.seed, dir) {
  signif(read_results(alpha, theta, n, p, data.seed, dir)$scalars$spi, 2)
}
