#' Bayesian Network Thresholding Heuristic Algorithm (BANTHA)
#'
#' This function provides a Bayesian network to summarize a network distribution
#' using the BANTHA method (Dahl, Bailey, Harris, Andros 2025). 
#'
#' @param samples One of three possible values:
#'   First option: An object of class \sQuote{array} containing posterior samples
#'   from a Bayesian network distribution. Each array element encodes one
#'   network as a binary adjacency matrix, where the rows and columns represent
#'   nodes. A value of 1 in cell \code{(i, j)} of an adjacency matrix indicates 
#'   a directed edge from node \code{i} to node \code{j}. The third dimension 
#'   of the array corresponds to different adjacency matrices (samples) from the 
#'   posterior distribution.  Second option: An object returned by the function
#'   \code{\link{gsh_loss_builder_new}} and used in calls to the function
#'   \code{\link{gsh_loss_builder_process}}.  Third option: An object returned by
#'   the function \code{\link{compute_expected_gsh_loss}}.
#' @param a A numeric scalar for the cost parameter of generalized Hamming
#'   distance used in the GSH loss.  The other cost parameter, \eqn{b}, is equal to
#'   \eqn{2 - a}.  Do not set unless \code{samples} above is the type from the first option.
#' @param n_candidates The number of possible starting states to use in calculating 
#'   the best matrix estimate. Using all samples as starting states is indicated
#'   by a value of zero.  Do not set unless \code{samples} above is the type from the first option.
#' @param candidates A three dimensional array (like that described for \code{samples})
#'   the best matrix estimate. Using all samples as starting states is indicated
#'   by a value of zero.  Do not set unless \code{samples} above is the type from the third option.
#' @param n_cores The number of CPU cores to use, i.e., the number of
#'   simultaneous calculations at any given time. A value of zero indicates to
#'   use all cores on the system.
#' @return A Bayes estimate in the form of an adjacency matrix found by minimizing 
#'   the Monte Carlo approximation of the expected GSH loss using the available 
#'   samples.
#'
#' @importFrom methods hasArg 
#' @export
#'
#' @references
#' D. B. Dahl, E. Bailey, D. Harris, J. Andros (2025), Bayesian Network Estimation from Posterior Samples,
#' <https://dahl.byu.edu/papers/dahl_bailey_harris_andros.pdf>.
#'
#' @examples
#' data(choi25)
#' bantha(choi25, a = 1)
#'
bantha <- function(samples, a = 1, n_candidates = 0, candidates = NULL, n_cores = 0) {
  if (is.array(samples)) {
    if (hasArg("candidates")) stop("Based on the type of the 'sample' argument, 'candidates' should not be provided.")
    .Call(.bantha, samples, a, n_candidates, n_cores)
  } else if (mode(samples) == "externalptr") {
    if (hasArg("a")) stop("Based on the type of the 'sample' argument, 'a' should not be provided.")
    if (hasArg("n_candidates")) stop("Based on the type of the 'sample' argument, 'n_candidates' should not be provided.")
    tag <- .Call(.get_external_ptr_tag, samples)
    if (tag == "gsh_builder") {
      if (hasArg("candidates")) stop("Based on the type of the 'sample' argument, 'candidates' should not be provided.")
      .Call(.bantha_big_data, samples, n_cores)
    } else if (tag == "gsh") {
      .Call(.bantha_psm, samples, candidates, n_cores)
    } else {
      stop("Unsupporrted external pointer type")
    }
  } else stop("Unsupported type for 'samples'")
}
