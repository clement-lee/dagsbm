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

#' Simulate from DAG-SBM
#'
#' @param alpha Real value less than 1, concentration parameter of the Pitman-Yor process
#' @param theta Real value, >= -alpha if alpha is between 0 (inclusive) and 1 (exclusive), or a positive integer multiple of |alpha| if alpha < 0
#' @param n Positive integer, number of nodes of the network i.e. order of corresponding graph
#' @param data.seed Positive integer, the seed for set.seed()
#' @param a0 Numeric, hyperparameter for prior of block matrix
#' @param b0 Numeric, hyperparameter for prior of block matrix
#' @param c0 Numeric, hyperparameter for prior of degree correction
#' @param d0 Numeric, hyperparameter for prior of degree correction
#' @importFrom igraph graph_from_adjacency_matrix as_adjacency_matrix
#' @importFrom stats rgamma rpois
#' @export
sim_dagsbm <- function(alpha, theta, n, data.seed, a0, b0, c0, d0) {
  ## sample Z from the Chinese restaurant process
  Z = rep(0,n)
  Z[1]=1
  set.seed(data.seed)
  for (i in 2:n){
    ## Compute frequencies in Z_{1:(i-1)}
    freq0 = as.numeric(table(Z[1:(i-1)]) )
    ## Labels of current blocks in Z_{1:(i-1)}
    labels = as.numeric(names(table(Z[1:(i-1)])))
    ## number of blocks in Z_{1:(i-1)}
    K = max(labels) ## or  K=max(Z[1:(i-1)])
    ## probabilities of i joining existing block or creating a new one
    probab = c(freq0-alpha,theta+alpha*K) 
    ## possible labels for i (old labels or K+1)
    poss_values = c(labels,K+1)
    ## Sample Z_i
    Z[i] = sample(poss_values,1,prob = probab)
  }
  ## Number of blocks in the simulated data
  K = max(Z)
  print(glue("Actual K: {K}"))
  ## sample C,  C_k,k' \sim Gamma(a,b)
  C= matrix(stats::rgamma(K^2,a0,b0), nrow = K)
  ## Degree corrected parameters, eta_i \sim Gamma(c,d)
  eta = rgamma(n,c0,d0)
  ## sample the topological order, sigma uniform over permutations of 1,...,n
  sigma = sample(1:n,n,replace=FALSE)
  phi = order(sigma)
  ## Adjacency matrix Y
  Y = matrix(rep(0,n^2),nrow=n)
  for (i in 1:n){
    for (j in 1:n){
      ## to include self-loop, set >= in the following condition
      ## if i is before j in the topological order, then there can be an edge
      if (phi[i]<phi[j]){
        ## Sample Y_{i,j} from Poisson
        Y[i,j] = stats::rpois(1,  eta[i]*eta[j]*C[Z[i],Z[j]])
      }
    }
  }
  ## return
  list(
    K = K,
    Y = Y |> igraph::graph_from_adjacency_matrix() |> igraph::as_adjacency_matrix(sparse = TRUE), # make sparse
    sigma = sigma - 1L # make 0-indexed
  )
}
