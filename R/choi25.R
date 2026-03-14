#' Posterior Samples from Choi Simulation
#'
#' Samples are provided from the simulated data used by Choi et al. (2020).
#'
#' @format An object of class \sQuote{list} containing 1,500 posterior samples
#'   from Choi et al. (2020)'s simulated posterior draws with 25 nodes. Each 
#'   list element contains a binary adjacency matrix, with nodes in the rows and
#'   columns indicating the presence or absence of directed edges.
#'
#' @usage data(choi25)
#'
#' @references
#' J. Choi, R. Chapkin, Y. Ni (2020). 
#' "Bayesian causal structural learning with zero-inflated poisson 
#' bayesian networks." Advances in Neural Information Processing Systems, 
#' 33: 5887–5897.
#'
"choi25"
