#' Sample the Number of Clusters from a Two-Parameter CRP
#'
#' Generates random samples of the number of clusters (tables) when sampling
#' from the two-parameter Chinese Restaurant Process (CRP), also known as the
#' Pitman-Yor process.
#'
#' @param n Integer specifying the number of samples to generate.
#' @param n_items Integer specifying the number of items (customers) in the set partition.
#' @param concentration Numeric value specifying the concentration (mass) parameter,
#'   must be greater than `-discount`. For the standard one-parameter CRP (Dirichlet
#'   process), this is typically denoted alpha and must be positive.
#' @param discount Numeric value in `[0, 1)` specifying the discount parameter.
#'   When `discount = 0`, this reduces to the standard one-parameter CRP.
#'
#' @return An integer vector of length `n` containing the sampled number of clusters.
#'
#' @details
#' The two-parameter CRP (Pitman-Yor process) generalizes the standard CRP by
#' introducing a discount parameter. The probability that a new customer sits at
#' a new table, given `k` existing tables and `i` customers already seated, is:
#'
#' \deqn{P(\text{new table}) = \frac{\alpha + k \cdot \sigma}{\alpha + i}}{P(new table) = (alpha + k * sigma) / (alpha + i)}
#'
#' where `alpha` is the concentration parameter and `sigma` is the discount parameter.
#'
#' The expected number of clusters for `n` items is approximately:
#' - For `discount = 0`: `alpha * log(n/alpha + 1)` (logarithmic growth)
#' - For `discount > 0`: `C * n^discount` for some constant `C` (polynomial growth)
#'
#' @seealso [dcrpk()] for the probability mass function.
#'
#' @example man/examples/rcrpk.R
#' @export
rcrpk <- function(n, n_items, concentration, discount) {
  .Call(.rcrpk, n, n_items, concentration, discount)
}
