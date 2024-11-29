#' Plot adjacency matrix with the rows and columns in topological order, using the edgelist and topological ordering
#'
#' @param df.edgelist Data frame of the edgelist of the graph
#' @param df.topo Data frame of the topological ordering of the vertices of the graph
#' @param weighted Boolean, is the graph edge-weighted (TRUE) or binary (FALSE)?
#' @importFrom dplyr select left_join
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot coord_cartesian theme_void theme element_rect geom_point aes scale_colour_gradient geom_abline
#' @keywords internal
plot_adj_mat <- function(df.edgelist, df.topo, weighted = TRUE) {
  weighted <- weighted && ("weight" %in% names(df.edgelist))
  if (weighted) {
    df0 <- df.edgelist |> dplyr::select(.data$from, .data$to, .data$weight)
  } else {
    df0 <- df.edgelist |> dplyr::select(.data$from, .data$to)
  }
  df0 <- df0 |>
  dplyr::left_join(df.topo, c("from" = "id")) |>
  dplyr::left_join(df.topo, c("to" = "id"))
  n <- nrow(df.topo)
  gg <- df0 |> 
    ggplot2::ggplot() +
    ggplot2::coord_cartesian(xlim = c(0, n), ylim = c(n, 0)) +
    ggplot2::theme_void(15) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(), 
      legend.position = "none",
      panel.background = ggplot2::element_rect(fill = "white")
    )
  if (weighted) { # for later use
    gg <- 
      gg + 
      ggplot2::geom_point(ggplot2::aes(.data$id_num.y, .data$id_num.x, col = .data$weight), size = 0.75) +
      ggplot2::scale_colour_gradient(low = "white", high = "black")
  } else {
    gg <- 
      gg + 
      ggplot2::geom_point(ggplot2::aes(.data$id_num.y, .data$id_num.x), size = 0.25, col = "black", alpha = 0.5)
  }
  gg +
    ggplot2::geom_abline(intercept = 0, slope = 1, col = "black", lty = "dashed")
}

#' Plot adjacency matrix with the rows and columns in topological order
#'
#' @param alpha Real value less than 1, concentration parameter of the Pitman-Yor process
#' @param theta Real value, >= -alpha if alpha is between 0 (inclusive) and 1 (exclusive), or a positive integer multiple of |alpha| if alpha < 0
#' @param n Positive integer, number of nodes of the network i.e. order of corresponding graph
#' @param p Real value between 0 and 1 inclusive, prior probability of finite regime
#' @param data.seed Positive integer, the seed for set.seed()
#' @param dir Boolean, is the SBM for directed graphs (TRUE) or directed acyclic graphs (FALSE)?
#' @param A_or_Y Boolean, is it the average (TRUE) or actual (FALSE) adjacency matrix that is required from the object of \code{\link{read_results}}?
#' @importFrom igraph graph_from_adjacency_matrix as_data_frame
#' @importFrom tibble as_tibble tibble
#' @importFrom glue glue
#' @importFrom ggplot2 labs theme element_text
#' @export
plot_adj_wrapper <- function(alpha, theta, n, p, data.seed, dir, A_or_Y) {
  alpha0 <- format(alpha, nsmall = 2)
  theta0 <- format(theta, nsmall = 2)
  p <- format(p, nsmall = 2)
  obj0 <- read_results(alpha0, theta0, n, p, data.seed, dir)
  df0.edgelist <- 
    (if (A_or_Y) obj0$A else obj0$Y) |> 
    igraph::graph_from_adjacency_matrix(weighted = "weight") |> 
    igraph::as_data_frame() |> tibble::as_tibble()
  df0.topo <- tibble::tibble(id = obj0$sigma.true+1, id_num = seq_along(.data$id))
  title0 <- 
    if (!A_or_Y) {
      glue::glue("alpha = {alpha}, theta = {theta}")
    } else {
      if (p == "0.00") {
        "infinite regime"
      } else if (p == "1.00") {
        "finite regime"
      } else {
        glue::glue("Pr(infinite) = {p}")
      }
    }
  plot_adj_mat(df0.edgelist, df0.topo, A_or_Y) +
    ggplot2::labs(title = title0) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 35))
}

#' Plot network with specific settings
#'
#' @param g igraph object
#' @param com Community memberships, of class "communities"
#' @param mar Scalar for the margins (all 4 sides of the bounding box)
#' @param seed Positive integer, random seed for the layout
#' @importFrom igraph layout_nicely
#' @export
plot_network <- function(g, com, mar, seed) {
  if (missing(com)) {
    set.seed(seed)
    plot(
      g,
      vertex.size = 7.5,
      vertex.label.cex = 0.01,
      edge.width = 0.25,
      edge.arrow.width = 0.75,
      edge.arrow.size = 0.25,
      margin = rep(mar, 4L),
      layout = igraph::layout_nicely
    )
  } else {
    set.seed(seed)
    plot(
      com, g,
      vertex.size = 7.5,
      vertex.label.cex = 0.01,
      edge.width = 0.25,
      edge.arrow.width = 0.75,
      edge.arrow.size = 0.25,
      margin = rep(mar, 4L),
      layout = igraph::layout_nicely
    )
  }
}

#' Plot posterior of topological ordering
#'
#' @param df Data frame to be passed to \code{\link{pivot_topo}}
#' @param colours Boolean, are colours to be added (default TRUE)?
#' @param mean.pos Boolean, are the mean positions to be included as points (default FALSE)?
#' @param low String, the colour for the lower end of the scale; ignored if colours is FALSE
#' @param true.pos Integer vector, a permutation of 1, 2, ..., length(true.pos), or NULL (default)
#' @importFrom ggplot2 ggplot labs scale_x_continuous scale_y_reverse theme_bw theme element_blank geom_tile aes scale_fill_gradient geom_point
#' @importFrom scales comma
#' @export
plot_topo <- function(df, colours = TRUE, mean.pos = FALSE, low = "yellow", true.pos = NULL) {
  ## df is NOT result of (but to be used in) pivot_topo()
  df1 <- df |> pivot_topo(v = true.pos)
  random <- (is.null(true.pos) || missing(true.pos))
  gg <- df1 |>
    ggplot2::ggplot() +
    ggplot2::labs(
      x = expression(paste("Position of vertex in ", sigma)),
      y = if (random) "An arbitrary topological ordering of vertices" else "True topological ordering"
    ) +
    ggplot2::scale_x_continuous(expand = c(0.0, 0.0)) +
    ggplot2::scale_y_reverse(expand = c(0.0, 0.0)) +
    ggplot2::theme_bw(22) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  if (colours) {
    gg <- gg +
      ggplot2::geom_tile(ggplot2::aes(.data$value, .data$position, fill = .data$density)) +
      ggplot2::scale_fill_gradient(
        low = low,
        high = "red",
        trans = "log",
        breaks = 2*10^(-1:-4),
        labels = scales::comma,
        )
  }
  if (random && mean.pos) {
    gg <- gg + ggplot2::geom_point(ggplot2::aes(.data$mean, .data$position), alpha = 0.1, size = 0.1)
  }
  gg
}

#' Wrapper of \code{\link{plot_topo}}
#'
#' @param alpha Real value less than 1, concentration parameter of the Pitman-Yor process
#' @param theta Real value, >= -alpha if alpha is between 0 (inclusive) and 1 (exclusive), or a positive integer multiple of |alpha| if alpha < 0
#' @param n Positive integer, number of nodes of the network i.e. order of corresponding graph
#' @param p Real value between 0 and 1 inclusive, prior probability of finite regime
#' @param data.seed Positive integer, the seed for set.seed()
#' @param dir Boolean, is the SBM for directed graphs (TRUE) or directed acyclic graphs (FALSE, default)?
#' @importFrom ggplot2 labs
#' @importFrom glue glue
#' @seealso \code{\link{obtain_spi}}, \code{\link{read_results}} and \code{\link{compute_bf}} with the same set of arguments
#' @export
plot_topo_wrapper <- function(alpha, theta, n, p, data.seed, dir = FALSE) {
  obj0 <- read_results(alpha, theta, n, p, data.seed, dir)
  plot_topo(obj0$phi, true.pos = order(obj0$sigma.true)) +
    ggplot2::labs(title = glue::glue("alpha = {alpha}, theta = {theta}"))
}

#' Plot adjacency matrix overlaid by similarity (co-clustering) matrix
#'
#' @param l List, usually returned by \code{\link{obtain_point_est}}
#' @param width Numeric, line width of the vertical & horizontal lines separating the blocks
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_reverse scale_fill_gradient geom_tile geom_point geom_hline geom_vline labs theme_void theme element_rect element_text unit
#' @export
plot_adj_psm <- function(l, width) {
  ## plot adj mat with co-clustering & clustering
  ## l is output from obtain_point_est()
  l$psm |>
  dplyr::filter(.data$Freq > 0) |>
  ggplot2::ggplot(ggplot2::aes(as.integer(.data$Var1), as.integer(.data$Var2))) +
    ggplot2::scale_x_continuous(expand = c(0.0, 0.0)) +
    ggplot2::scale_y_reverse(expand = c(0.0, 0.0)) +
    ggplot2::scale_fill_gradient(low = "white", high = "red", name = "co-clustering probability") +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$Freq), stat = "identity") +
    ggplot2::geom_point(data = l$adj |> dplyr::filter(.data$Freq > 0), size = 0.1, alpha = 0.5) +
    ggplot2::geom_hline(yintercept = c(0, l$cN) + 0.5, col = 4, lty = 2, lwd = width, alpha = 0.5) +
    ggplot2::geom_vline(xintercept = c(0, l$cN) + 0.5, col = 4, lty = 2, lwd = width, alpha = 0.5) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void(15) +
    ggplot2::theme(panel.border = ggplot2::element_rect(),
                   legend.position = "top",
                   legend.title = ggplot2::element_text(size = 20),
                   legend.text = ggplot2::element_text(angle = 45, size = 20, vjust = 0.45),
                   legend.key.size = ggplot2::unit(1, "cm"))
}

#' Plot subset of the simulation results
#'
#' @param df Data frame of MCMC samples
#' @param alpha Real value less than 1, concentration parameter of the Pitman-Yor process
#' @param theta Real value, >= -alpha if alpha is between 0 (inclusive) and 1 (exclusive), or a positive integer multiple of |alpha| if alpha < 0
#' @param data.seeds Integer vector, the seeds used for set.seed()
#' @param ns Integer vector, numbers of nodes of the network i.e. order of corresponding graph; must be of same length as data.seeds
#' @param scales String that is accepted by the argument of the same name in ggplot2::facet_grid
#' @importFrom dplyr filter semi_join mutate vars
#' @importFrom rlang .data
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes geom_vline facet_grid labs scale_x_continuous scale_y_continuous theme_bw theme element_text
#' @importFrom ggpattern geom_bar_pattern scale_pattern_manual
#' @importFrom glue glue
#' @export
plot_subset <- function(df, alpha, theta, data.seeds, ns = c(250, 500, 1000), scales = "free") {
  if (length(data.seeds) != length(ns)) {
    stop("plot_subset: lengths of data.seeds & ns have to be equal.")
  }
  df |>
  dplyr::filter(
    .data$alpha.true == {{ alpha }}, 
    .data$theta.true == {{ theta }}, 
    .data$p_fin %in% c(0.0, 1.0)
  ) |>
  dplyr::semi_join(
    tibble::tibble(n = ns, data.seed = data.seeds),
    by = c("n", "data.seed")
  ) |>
  dplyr::mutate(regime = ifelse(.data$p_fin == 0.0, "infinite", "finite")) |>
  ggplot2::ggplot() +
    ggpattern::geom_bar_pattern(ggplot2::aes(.data$K, fill = .data$regime, pattern = .data$regime), position = "dodge") +
    ggpattern::scale_pattern_manual(values = c(infinite = "stripe", finite = "none")) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = .data$K.true), col = 2, lty = 2) +
    ggplot2::facet_grid(dplyr::vars(.data$model), dplyr::vars(.data$nfull), scales = scales) +
    ggplot2::labs(title = glue::glue("alpha = {alpha}, theta = {theta}"), x = NULL) +
    ggplot2::scale_x_continuous(minor_breaks = NULL) +
    ggplot2::scale_y_continuous(minor_breaks = NULL) +
    ggplot2::theme_bw(12) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12))
}

#' Trace plots of MCMC
#'
#' @param df Data frame of MCMC samples
#' @importFrom ggplot2 ggplot geom_line aes facet_wrap labs theme_bw
#' @importFrom dplyr group_by mutate
#' @importFrom rlang .data
#' @export
plot_trace <- function(df) {
  n.var <- df$name |> unique() |> length()
  df |>
  dplyr::group_by(.data$name) |>
  dplyr::mutate(index = seq_along(.data$value)) |>
  ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(.data$index, .data$value)) +
    ggplot2::facet_wrap(~name, n.var, 1, "free") +
    ggplot2::labs(y = NULL, x = "Iteration") +
    ggplot2::theme_bw(25)
}

#' Plot posterior density
#'
#' @param df Data frame of MCMC samples
#' @importFrom ggplot2 ggplot geom_density aes facet_wrap labs theme_bw
#' @importFrom rlang .data
#' @export
plot_density <- function(df) {
  n.var <- df$name |> unique() |> length()
  df |>
  ggplot2::ggplot() +
    ggplot2::geom_density(ggplot2::aes(.data$value)) +
    ggplot2::facet_wrap(~name, n.var, 1, "free") +
    ggplot2::labs(x = "value", y = NULL) +
    ggplot2::theme_bw(25)
}
