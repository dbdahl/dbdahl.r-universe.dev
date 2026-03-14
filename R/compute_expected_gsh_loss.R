#' Compute Expected Generalized Structural Hamming (GSH) Loss
#'
#' This function computes the expected GSH loss (Dahl, Bailey, Harris, Andros, 2025)
#' for a given adjacency matrix.
#'
#' @param network A Bayesian network represented by an adjacency matrix.
#' @param gsh A pointer to an initialized GSH structure containing the posterior 
#'   probabilities of edges between each pair of nodes, which can be found using
#'   the \code{\link{initialize_expected_gsh_loss}} function.
#' @return The expected GSH loss for the given Bayesian network based on the posterior 
#'   probabilities of edges.
#'
#' @export
#'
#' @references
#' D. B. Dahl, E. Bailey, D. Harris, J. Andros (2025), Bayesian Network Estimation from Posterior Samples,
#' <https://dahl.byu.edu/papers/dahl_bailey_harris_andros.pdf>.
#'
#' @examples
#' data(choi25)
#' gsh <- initialize_expected_gsh_loss(choi25, a = 1)
#' compute_expected_gsh_loss(choi25[,,1], gsh)
#'
compute_expected_gsh_loss <- function(network, gsh) {
  .Call(.compute_expected_gsh_loss, network, gsh)
}
