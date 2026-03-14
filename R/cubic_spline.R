#' Create a Natural Cubic Spline
#'
#' Constructs a natural cubic spline for smooth interpolation on the interval \code{[0,1]}.
#' The spline passes through all control points and has zero second derivative
#' at the endpoints (natural boundary conditions). Cubic splines can be passed as
#' the `concentration`, `log_skew`, or `tail_shape` arguments to [rlorenzip()] and
#' [dlorenzip()] so that these parameters vary smoothly across sequential steps.
#'
#' @param x Numeric vector of control point x-coordinates. Must be strictly increasing
#'   with `x[1] = 0` and `x[length(x)] = 1`.
#' @param y Numeric vector of control point y-values. Must have the same length as `x`.
#'
#' @return An external pointer object of class `"cubic_spline"` representing
#'   the cubic spline curve.
#'
#' @seealso [spline_evaluate()] for evaluating the spline at specific points.
#'
#' @examples
#' # Linear function (just endpoints)
#' spline <- cubic_spline(c(0, 1), c(0, 1))
#' spline_evaluate(spline, 0.5) # Returns 0.5
#'
#' # Curved function with intermediate points
#' spline <- cubic_spline(c(0, 0.25, 0.75, 1), c(0, 0.1, 0.9, 1))
#' spline_evaluate(spline, seq(0, 1, by = 0.1))
#'
#' @export
cubic_spline <- function(x, y) {
  .Call(.new_cubic_spline_extptr, x, y)
}

#' Evaluate a Cubic Spline
#'
#' Evaluates a cubic spline at specified points in \code{[0, 1]}.
#'
#' @param spline A cubic spline object created by [cubic_spline()].
#' @param x Numeric vector of points in \code{[0, 1]} at which to evaluate the spline.
#'
#' @return Numeric vector of the same length as `x` containing the spline
#'   values at the specified points.
#'
#' @seealso [cubic_spline()] for creating cubic spline objects.
#'
#' @examples
#' # Create a spline
#' spline <- cubic_spline(c(0, 0.5, 1), c(0, 0.8, 1))
#'
#' # Evaluate at specific points
#' spline_evaluate(spline, c(0, 0.25, 0.5, 0.75, 1))
#'
#' @export
spline_evaluate <- function(spline, x) {
  UseMethod("spline_evaluate")
}

#' @export
spline_evaluate.cubic_spline <- function(spline, x) {
  .Call(.cubic_spline_evaluate, spline, x)
}

#' @export
spline_evaluate.default <- function(spline, x) {
  stop("spline_evaluate() requires a cubic_spline object")
}

#' @export
print.cubic_spline <- function(x, ...) {
  .Call(.print_cubic_spline, x)
  invisible(x)
}

#' Plot a Cubic Spline
#'
#' Plots a cubic spline curve.
#'
#' @param x A cubic spline object (functions as the `spline` parameter).
#' @param n_points Number of points to use for plotting (default 201).
#' @param ... Additional arguments passed to [plot()].
#'
#' @return Invisibly returns the cubic spline object.
#'
#' @examples
#' # Create and plot a curved spline
#' spline <- cubic_spline(c(0, 0.2, 0.5, 0.8, 1), c(5, 10, 25, 40, 50))
#' plot(spline, main = "Concentration Curve")
#'
#' @export
plot.cubic_spline <- function(x, n_points = 201L, ...) {
  t <- seq(0, 1, length.out = n_points)
  y <- spline_evaluate(x, t)
  plot(
    t,
    y,
    type = "l",
    xlab = "Position (0 to 1)",
    ylab = "Value",
    ...
  )
  invisible(x)
}
