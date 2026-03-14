n_items <- 20
weights <- c(1, 2, 3, 4)
concentration <- 5
partition <- c(2, 4, 6, 8)

diag <- lorenzip_diagnostics(
  partition = partition,
  n_items = n_items,
  target = weights,
  concentration = concentration
)

diag$summary
diag$steps
