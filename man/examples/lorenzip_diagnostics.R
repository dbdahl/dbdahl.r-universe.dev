n_items <- 20
weights <- c(1, 2, 3, 4)
concentration <- 5
partition <- c(2, 4, 6, 8)

diag_path <- lorenzip_diagnostics(
  partition = partition,
  n_items = n_items,
  target = weights,
  concentration = concentration
)

diag_path$summary
diag_path$steps

diag_mc <- lorenzip_diagnostics(
  n_samples = 50,
  n_items = n_items,
  target = weights,
  concentration = concentration
)

diag_mc$summary
head(diag_mc$steps)
