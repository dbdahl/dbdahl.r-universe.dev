#' Plan Splinters for the Splintered Clustering Procedure
#'
#' This function provides a randomize splintering plan for the splintered
#' clustering procedure. Note that \sQuote{nSplinters} * \sQuote{splinterSize}
#' must be at least \sQuote{nItems}.
#'
#' @param n_items Total number of items to cluster.
#' @param n_splinters Number of splinters to form.
#' @param splinter_size Number of items in each splinter.
#' @param n_anchors Number of anchors, i.e., items found in each splinter.
#'
#' @return A list containing \sQuote{nSplinters} lists, each containing a vector
#'   \sQuote{indices} giving the item numbers for a given splinter.
#' @export
#'
#' @examples
#' plan(20, 5, 8)
#'
plan <- function(
  n_items,
  n_splinters,
  splinter_size,
  n_anchors = 0
) {
  if (n_anchors == 0) {
    plan_engine(n_items, n_splinters, splinter_size, identity)
  } else {
    anchors <- sample.int(n_items, n_anchors)
    nonAnchors <- setdiff(seq_len(n_items), anchors)
    mapper <- function(x) c(anchors, nonAnchors[x])
    plan_engine(
      n_items - n_anchors,
      splinter_size - n_anchors,
      n_splinters,
      mapper
    )
  }
  #
}

plan_engine <- function(
  n_items,
  n_splinters,
  splinter_size,
  mapper = identity
) {
  results <- list()
  index <- 1
  counts <- rep(0, n_items)
  for (index in seq_len(n_splinters)) {
    p <- 1 * (counts == min(counts))
    x <- if (sum(p) < splinter_size) {
      equal_min <- counts == min(counts)
      p <- 1 * (!equal_min)
      c(
        which(equal_min),
        sample.int(n_items, splinter_size - sum(equal_min), prob = p)
      )
    } else {
      sample.int(n_items, splinter_size, prob = p)
    }
    counts[x] <- counts[x] + 1
    results[[index]] <- list(indices = mapper(x))
    index <- index + 1
  }
  results
}
