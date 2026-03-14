#' Create a Piecewise Linear Lorenz Curve
#'
#' Constructs a piecewise linear Lorenz curve from an empirical partition.
#' The curve exactly interpolates the cumulative distribution of the partition.
#'
#' @inheritParams params
#'
#' @return An external pointer object of class `c("lorenz_linear", "lorenz")`
#'   representing the Lorenz curve.
#'
#' @seealso [lorenz_ispline()] for smooth I-spline Lorenz curves,
#'   [gini()] for computing the Gini coefficient,
#'   [lorenz_increments()] for extracting Lorenz curve increments.
#'
#' @example man/examples/lorenz.R
#' @export
lorenz_linear <- function(partition) {
  .Call(.new_lorenz_linear_extptr, partition)
}

#' Create an I-Spline Lorenz Curve
#'
#' Constructs a smooth I-spline Lorenz curve fitted to an empirical partition.
#' The curve is constrained to be monotonically increasing and pass through
#' (0, 0) and (1, 1).
#'
#' @param n_interior_knots Number of interior knots for the I-spline basis.
#'   If `NULL` (default), automatically set to `max(1, min(length(partition) - 2, 10))`.
#' @param degree Degree of the I-spline (default 2, giving cubic splines).
#' @param lambda Roughness penalty weight (default 1e-6). Larger values give
#'   smoother curves.
#'
#' @inheritParams params
#'
#' @return An external pointer object of class `c("lorenz_ispline", "lorenz")`
#'   representing the Lorenz curve.
#'
#' @seealso [lorenz_linear()] for piecewise linear Lorenz curves,
#'   [gini()] for computing the Gini coefficient,
#'   [lorenz_increments()] for extracting Lorenz curve increments.
#'
#' @example man/examples/lorenz.R
#' @export
lorenz_ispline <- function(
  partition,
  n_interior_knots = NULL,
  degree = 2L,
  lambda = 1e-6
) {
  .Call(
    .new_lorenz_ispline_extptr,
    partition,
    n_interior_knots,
    degree,
    lambda
  )
}

#' Evaluate a Lorenz Curve
#'
#' Evaluates a Lorenz curve at specified points in \[0, 1\].
#'
#' @param x Numeric vector of points in \[0, 1\] at which to evaluate the curve.
#'
#' @inheritParams params
#'
#' @return Numeric vector of the same length as `x` containing the Lorenz curve
#'   values L(x), representing the cumulative share of mass held by the poorest
#'   x fraction.
#'
#' @example man/examples/lorenz.R
#' @export
lorenz_evaluate <- function(lorenz, x) {
  UseMethod("lorenz_evaluate")
}

#' @export
lorenz_evaluate.lorenz_linear <- function(lorenz, x) {
  .Call(.lorenz_evaluate_linear, lorenz, x)
}

#' @export
lorenz_evaluate.lorenz_ispline <- function(lorenz, x) {
  .Call(.lorenz_evaluate_ispline, lorenz, x)
}

#' Compute the Gini Coefficient
#'
#' Computes the Gini coefficient from a Lorenz curve or directly from a partition.
#' The Gini coefficient measures inequality, ranging from 0 (perfect equality)
#' to 1 (maximum inequality).
#'
#' @param x Either a Lorenz curve object created by [lorenz_linear()] or
#'   [lorenz_ispline()], or a numeric vector representing a partition.
#'
#' @return A single numeric value representing the Gini coefficient.
#'
#' @example man/examples/lorenz.R
#' @export
gini <- function(x) {
  UseMethod("gini")
}

#' @export
gini.lorenz_linear <- function(x) {
  .Call(.lorenz_gini_linear, x)
}

#' @export
gini.lorenz_ispline <- function(x) {
  .Call(.lorenz_gini_ispline, x)
}

#' @export
gini.numeric <- function(x) {
  .Call(.gini_coefficient, x)
}

#' @export
gini.default <- function(x) {
  stop("gini() requires a lorenz object or numeric partition vector")
}

#' Compute the Gini Coefficient from a Partition
#'
#' Computes the Gini coefficient directly from a partition vector without
#' creating a Lorenz curve object.
#'
#' @inheritParams params
#'
#' @return A single numeric value representing the Gini coefficient.
#'
#' @seealso [gini()] for the generic function that also works with Lorenz curve
#'   objects.
#'
#' @example man/examples/lorenz.R
#' @export
gini_coefficient <- function(partition) {
  .Call(.gini_coefficient, partition)
}

#' Extract Increments from a Lorenz Curve
#'
#' Extracts the Lorenz curve increments
#' \eqn{\omega_j = L(j/k) - L((j-1)/k)} for a specified number of clusters
#' \eqn{k}. These are the target mean profile weights used by the Lorenz-IP
#' distribution, in non-decreasing order (smallest to largest).
#'
#' @param n_clusters Positive integer specifying the number of clusters.
#'
#' @inheritParams params
#'
#' @return Numeric vector of length `n_clusters` summing to 1, in
#'   non-decreasing order.
#'
#' @example man/examples/lorenz.R
#' @export
lorenz_increments <- function(lorenz, n_clusters) {
  UseMethod("lorenz_increments")
}

#' @export
lorenz_increments.lorenz_linear <- function(lorenz, n_clusters) {
  .Call(.lorenz_weights_linear, lorenz, n_clusters)
}

#' @export
lorenz_increments.lorenz_ispline <- function(lorenz, n_clusters) {
  .Call(.lorenz_weights_ispline, lorenz, n_clusters)
}

#' @export
print.lorenz_linear <- function(x, ...) {
  .Call(.print_lorenz_linear, x)
  invisible(x)
}

#' @export
print.lorenz_ispline <- function(x, ...) {
  .Call(.print_lorenz_ispline, x)
  invisible(x)
}

#' Plot a Lorenz Curve
#'
#' Plots a Lorenz curve with the line of equality (45-degree line) for reference.
#'
#' @param x A Lorenz curve object (functions as the `lorenz` parameter).
#' @param n_points Number of points to use for plotting (default 101).
#' @param ... Additional arguments passed to [plot()].
#'
#' @return Invisibly returns the Lorenz curve object.
#'
#' @example man/examples/lorenz.R
#' @importFrom graphics abline
#' @export
plot.lorenz <- function(x, n_points = 101L, ...) {
  p <- seq(0, 1, length.out = n_points)
  l <- lorenz_evaluate(x, p)
  plot(
    p,
    l,
    type = "l",
    xlab = "Cumulative share of clusters",
    ylab = "Cumulative share of items",
    xlim = c(0, 1),
    ylim = c(0, 1),
    ...
  )
  abline(0, 1, lty = 2, col = "gray")
  invisible(x)
}
