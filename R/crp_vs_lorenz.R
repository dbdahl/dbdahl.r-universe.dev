#' Monte Carlo CRP vs Lorenz Diagnostics
#'
#' Approximates total variation and KL divergence between the two-parameter CRP
#' and the Lorenz partition distribution on set partitions using Monte Carlo samples from the CRP.
#' The KL direction is CRP relative to Lorenz.
#'
#' @param crp_concentration Numeric value specifying the CRP concentration parameter.
#'   Must be greater than `-crp_discount`.
#' @param crp_discount Numeric value in `[0, 1)` specifying the CRP discount parameter.
#' @param n_samples Integer specifying the number of Monte Carlo samples used in
#'   the approximation.
#' @param k_weights Optional numeric vector of non-negative weights for the prior
#'   on `k` used by the Lorenz distribution. Entry `j` is proportional to Pr(`k = j`).
#'   Weights are renormalized after truncation to `k <= n_items`. If `NULL`, the
#'   Lorenz prior on `k` is taken to match the CRP implied distribution (unless the
#'   Lorenz parameters fix `k`).
#' @param log_k_weights Logical; if `TRUE`, interpret `k_weights` as log weights
#'   or log probabilities.
#' @inheritParams params
#'
#' @return A list with elements `tv` and `kl` containing the Monte Carlo estimates.
#'
#' @details
#' The comparison is performed on integer partitions (size profiles), which is
#' equivalent to comparing the induced exchangeable set partition distributions.
#' If the Lorenz distribution assigns zero probability to any integer partition with
#' positive CRP mass, the KL estimate will be infinite.
#'
#' @example man/examples/crp_vs_lorenz.R
#' @export
crp_vs_lorenz <- function(
  n_items,
  target,
  concentration,
  log_skew = 0,
  tail_shape = NULL,
  buffer_c = 1,
  crp_concentration,
  crp_discount,
  n_samples = 1000L,
  k_weights = NULL,
  log_k_weights = FALSE
) {
  .Call(
    .crp_vs_lorenz,
    n_items,
    target,
    concentration,
    log_skew,
    tail_shape,
    buffer_c,
    crp_concentration,
    crp_discount,
    n_samples,
    k_weights,
    log_k_weights
  )
}
