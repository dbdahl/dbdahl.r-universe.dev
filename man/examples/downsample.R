n_items <- 100
k_weights <- c(0, 1, 1) # Equal probability that k=2 and k=3
target <- lorenz_ispline(c(1, 2, 4))
concentration <- 1

downsample(
  n_items,
  k_weights,
  target,
  concentration = concentration,
  n_samples = 200
)
