#' Add Bayesian Networks to GSH Builder Structure
#'
#' This function adds a Bayesian network represented by an adjacency matrix to an existing GSH 
#' builder structure.
#'
#' @param x A network as an adjacency matrix.
#' @param gsh_loss_builder A pointer to a GSH builder structure.
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
#' gsh_loss_builder_process(choi25[,,1], gsh_builder)
#'
gsh_loss_builder_process <- function(x, gsh_loss_builder) {
  .Call(.gsh_loss_builder_process, x, gsh_loss_builder)
}
