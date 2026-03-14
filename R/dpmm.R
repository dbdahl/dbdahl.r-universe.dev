#' Generate Data from the Neal (2000) Model
#'
#' This function generates data from the model in Section 8 of "Markov Chain Sampling Methods for Dirichlet Process Mixture Models" by Radford M. Neal,
#' Journal of Computational and Graphical Statistics, Vol. 9, No. 2. (Jun., 2000), pp. 249-265.
#'
#' @return A list giving the true partition, the associated cluster parameters, and the simulated data.
#'
#' @export
#' @examples
#' dpmm_generate(10)
#'
dpmm_generate <- function(n_items, mass = 1.0, sd = 0.1) {
  .Call(.dpmm_generate, n_items, mass, sd)
}

#' Sample from the Posterior Partition Distribution of the Neal (2000) Model
#'
#' This function samples from the posterior partition distribution of the Neal (2000) Model.  See Section 8 of
#' "Markov Chain Sampling Methods for Dirichlet Process Mixture Models" by Radford M. Neal,
#' Journal of Computational and Graphical Statistics, Vol. 9, No. 2. (Jun., 2000), pp. 249-265.
#'
#' @return A matrix in which each column represents the partition of the items.
#'
#' @export
#' @examples
#' x <- dpmm_generate(10)
#' dpmm_fit(x$data, 100, 10, 2)
#'
dpmm_fit <- function(data, n_samples, burnin, thin, mass = 1.0, sd = 0.1, algorithm = "neal8") {
  if (is.list(data)) {
    .Call(.dpmm_fit3, data, n_samples, burnin, thin, mass, sd)
  } else if (algorithm == "neal3") {
    .Call(.dpmm_fit, data, n_samples, burnin, thin, mass, sd)
  } else if (algorithm == "neal8") {
    .Call(.dpmm_fit2, data, n_samples, burnin, thin, mass, sd)
  } else {
    stop("Unrecognized algorithm")
  }
}

#' @export
neal_data <- function() {
    .Call(.neal_data)
}
