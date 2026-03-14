#' Probability Mass Function of the Number of Clusters from a Two-Parameter CRP
#'
#' Evaluates the probability mass function (PMF) for the number of clusters (tables)
#' when sampling from the two-parameter Chinese Restaurant Process (CRP), also known
#' as the Pitman-Yor process.
#'
#' @param x Integer vector of values at which to evaluate the PMF. Values must be
#'   between 1 and `n_items` inclusive; values outside this range return 0 (or `-Inf`
#'   if `log = TRUE`).
#' @param n_items Integer specifying the number of items (customers) in the set partition.
#' @param concentration Numeric value specifying the CRP concentration (mass)
#'   parameter alpha, must be greater than `-discount`. For the standard
#'   one-parameter CRP (Dirichlet process), alpha must be positive. Note: this
#'   is the CRP parameter alpha, not the Lorenz-IP concentration gamma.
#' @param discount Numeric value in `[0, 1)` specifying the discount parameter.
#'   When `discount = 0`, this reduces to the standard one-parameter CRP.
#' @param log Logical; if `TRUE`, returns log-probabilities. Default is `FALSE`.
#'
#' @return A numeric vector of the same length as `x` containing the (log-)probabilities.
#'
#' @details
#' The PMF for the number of clusters `K` in a sample of `n` items from a
#' two-parameter CRP with concentration `alpha` and discount `sigma` is:
#'
#' \deqn{P(K = k | n, \alpha, \sigma) = \frac{S_\sigma(n, k) \cdot C(\alpha, \sigma, k)}{(\alpha+1)_{n-1}}}{P(K = k | n, alpha, sigma) = S_sigma(n, k) * C(alpha, sigma, k) / (alpha+1)_(n-1)}
#'
#' where:
#' - `S_sigma(n, k)` is the generalized unsigned Stirling number of the first kind
#' - `C(alpha, sigma, k) = (alpha + sigma) * (alpha + 2*sigma) * ... * (alpha + (k-1)*sigma)`
#' - `(alpha+1)_(n-1)` is the rising factorial (Pochhammer symbol)
#'
#' For `discount = 0` (standard CRP), this simplifies to:
#' \deqn{P(K = k | n, \alpha) = \frac{|s(n, k)| \cdot \alpha^{k-1}}{(\alpha+1)_{n-1}}}{P(K = k | n, alpha) = |s(n, k)| * alpha^(k-1) / (alpha+1)_(n-1)}
#'
#' where `|s(n, k)|` is the unsigned Stirling number of the first kind.
#'
#' All computations are performed in log-space for numerical stability, making
#' this function suitable for large `n_items` (tested up to 100,000).
#'
#' @seealso [rcrpk()] for sampling from this distribution.
#'
#' @example man/examples/dcrpk.R
#' @export
dcrpk <- function(x, n_items, concentration, discount, log = FALSE) {
  .Call(.dcrpk, x, n_items, concentration, discount, log)
}
