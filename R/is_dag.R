#' Check if a Graph is a Directed Acyclic Graph (DAG)
#'
#' This function checks whether the provided graph is a Directed Acyclic Graph (DAG).
#
#' @param candidate An adjacency matrix representing a candidate graph.
#' @return A logical value indicating whether or not the graph is a DAG.
#'
#' @export
#'
#' @examples
#' # Example 1: This demonstration checks if a real DAG returns true.
#' dag_matrix <- matrix(c(0, 1, 0,
#'                        0, 0, 1,
#'                        0, 0, 0), nrow = 3, byrow = TRUE)
#' is_dag(dag_matrix)
#'
#' # Example 2: This demonstration checks if a graph with a cycle returns false.
#' cycle_matrix <- matrix(c(0, 1, 0,
#'                          0, 0, 1,
#'                          1, 0, 0), nrow = 3, byrow = TRUE)
#' is_dag(cycle_matrix)
#'
is_dag <- function(candidate) {
  .Call(.is_dag, candidate)
}
