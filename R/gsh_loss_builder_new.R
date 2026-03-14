#' Initialize Expected Generalized Structural Hamming (GSH) Loss Builder
#'  
#' Initialize an expected GSH Loss builder to process samples one at a time.
#'
#' @param n_items Number of nodes in each network to process.
#' @param n_samples Number of networks to process.
#' @param n_candidates Number of networks to sweeten.
#' @param a A numeric scalar for the cost parameter of generalized Hamming
#'   distance used in the GSH loss.  The other cost parameter, \eqn{b}, is equal to
#'   \eqn{2 - a}.
#' @return A pointer to an initialized GSH builder structure to process samples.
#'
#' @export
#'
#' @references
#' D. B. Dahl, E. Bailey, D. Harris, J. Andros (2025), Bayesian Network Estimation from Posterior Samples,
#' <https://dahl.byu.edu/papers/dahl_bailey_harris_andros.pdf>.
#'
#' @examples
#' data(choi25)
#' gsh_builder <- gsh_loss_builder_new(ncol(choi25), nrow(choi25), 4, a = 1)
#'
gsh_loss_builder_new <- function(n_items, n_samples, n_candidates, a = 1) {
  .Call(.gsh_loss_builder_new, n_items, n_samples, n_candidates, a)
}
