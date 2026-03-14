#' @export
new_tree_beta_gaussian <- function(depth, alpha, beta, mean, sd) {
  .Call(.new_tree_beta_gaussian, depth, alpha, beta, mean, sd)
}

#' @export
new_tree_from_funcs <- function(depth, break_fn, atom_fn) {
  .Call(.new_tree_from_funcs, depth, break_fn, atom_fn)
}

#' @export
new_dm_beta_gaussian <- function(depth, alpha, beta, mean, sd, method = "direct") {
  .Call(.new_dm_beta_gaussian, depth, alpha, beta, mean, sd, method)
}

#' @export
new_dm_from_breaks_and_atoms <- function(breaks, atoms) {
  .Call(.new_dm_from_breaks_and_atoms, breaks, atoms)
}

#' @export
nonzero_leaves <- function(tree) {
  .Call(.nonzero_leaves, tree)
}

#' @export
print.treesb <- function(tree, include_zeros = FALSE) {
  invisible(.Call(.print_treesb, tree, include_zeros))
}

#' @export
to_discrete_measure <- function(ptr) {
  .Call(.to_discrete_measure, ptr)
}

#' @export
dm_sample_partition <- function(ptr, n_items, n_samples) {
  .Call(.dm_sample_partition, ptr, n_items, n_samples)
}

#' @export
dm_sample <- function(ptr, n_samples) {
  .Call(.dm_sample, ptr, n_samples)
}

#' @export
dm_atoms <- function(ptr) {
  .Call(.dm_atoms, ptr)
}

#' @export
dm_weights <- function(ptr) {
  .Call(.dm_weights, ptr)
}

#' @export
path_to_leaf <- function(leaf_index, depth) {
  .Call(.path_to_leaf_r, leaf_index, depth)
}

#' @export
path_to_internal <- function(internal_index, depth) {
  .Call(.path_to_internal_r, internal_index, depth)
}


#' @export
dm_weights_from_breaks <- function(breaks) {
  .Call(.dm_weights_from_breaks_r, breaks)
}
