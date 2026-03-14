#' Common Parameters for lorenzpartition
#'
#' @param n_items Integer specifying the total number of items to partition.
#' @param target Either a numeric vector of positive weights
#'   `w_1, ..., w_k` proportional to the target mean profile
#'   `omega_j = L(j/k) - L((j-1)/k)` of the paper (the Lorenz-curve
#'   increments, not the cumulative Lorenz values `L(j/k)`).
#'   **Weights are automatically sorted to non-decreasing order** and
#'   **normalized internally to probabilities**; any positive scaling is
#'   allowed. Or, a Lorenz curve object created by [lorenz_linear()] or
#'   [lorenz_ispline()]. When a numeric vector is provided, `n_clusters`
#'   (number of clusters) is determined by its length. Use
#'   [lorenz_increments()] to extract the increments explicitly.
#' @param n_clusters Integer specifying the number of clusters. Required when
#'   `k` cannot be inferred from the lengths of `target`, `concentration`,
#'   `log_skew`, or `tail_shape`. Typically needed when `target` is a Lorenz
#'   curve object and all other parameters are scalars.
#'   When `NULL`, behavior depends on the function (see individual function documentation).
#' @param concentration Either a positive scalar, a numeric vector of length k-1
#'   (one per sequential step), or a cubic spline object created by [cubic_spline()]
#'   specifying the concentration parameters (`gamma_j` in the paper). Higher values concentrate the cluster
#'   sizes closer to the target weights times the number of items. When a scalar,
#'   the same value is used for every step. When a cubic spline is provided, it is
#'   evaluated at k-1 equally-spaced points (1/k, 2/k, ..., (k-1)/k).
#' @param log_skew Log-skew parameter (eta in the paper). Can be a numeric scalar,
#'   a numeric vector of length k-1, or a cubic spline object created by
#'   [cubic_spline()]. Missing or `NULL` defaults to 0. When a cubic spline is
#'   provided, it is evaluated at k-1 equally-spaced points (1/k, 2/k, ..., (k-1)/k).
#'   Controls left/right skew in the TiDaL or TaDPoLe kernel.
#' @param tail_shape Tail-shape parameter (alpha in the paper) for TaDPoLe. Can be
#'   a numeric scalar, a numeric vector of length k-1, or a cubic spline object
#'   created by [cubic_spline()]. If `NULL` (default), the TiDaL kernel is used.
#'   When provided, TaDPoLe is used and `tail_shape` controls tail heaviness.
#' @param buffer_c Positive scalar buffer scale `c` for the mean-space band.
#'   The per-step buffer is `delta_j = c / (1 + gamma_j)` and is capped at half
#'   the available support width. Larger values enforce more distance from the
#'   support boundaries. The default `buffer_c = 1` matches the paper's
#'   definition `delta_j = 1 / (1 + gamma_j)`.
#' @param log If `TRUE`, returns the natural logarithm. If `FALSE` (default),
#'   returns the value directly.
#' @param partition Numeric vector of positive values representing cluster sizes
#'   or weights. The partition does not need to be sorted; functions that take
#'   partitions sort to non-decreasing order internally.
#' @param lorenz A Lorenz curve object created by [lorenz_linear()] or
#'   [lorenz_ispline()].
#'
#' @name params
#' @keywords internal
NULL
