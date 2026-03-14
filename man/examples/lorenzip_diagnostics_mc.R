n_items <- 1000
weights <- c(2, 5, 10, 25, 50)
plot(lorenz_ispline(weights))

concentration <- c(9, 6, 5, 1)
log_skew <- log(1.2 / 0.8)

diag <- lorenzip_diagnostics(
  n_samples = 20000,
  n_items = n_items,
  target = weights,
  concentration = concentration,
  log_skew = log_skew
)

diag$summary
diag$steps
