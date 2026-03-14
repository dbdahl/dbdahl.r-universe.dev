#' Computed Expected Loss from Splinters
#'
#' This function computes the expected loss from splinters.
#'
#' @inheritParams splinclust
#' @param estimate A numeric vector giving an estimated clustering.
#'
#' @return A numeric vector.
#'
#' @export
#' @examples
#' cat("To do.\n")
#'
splinclust_expected_loss <- function(estimate, splinters, useVI=TRUE, a=1.0) {
  out <- check_splinters(splinters)
  len2 <- length(estimate)
  if (len2 != out$nItems) stop(sprintf("Length of estimate (%s) does not match that implied by splinters (%s).", len2, nItems))
  if (a < 0 || a > 2) stop("'a' should be in [0,2].")
  .Call(.expected_loss, estimate, out$ptrs, useVI, a)
}
