n_items <- 30
k <- 4
curve <- lorenz_ispline(c(1, 5, 15, 45))
plot(curve)
concentration <- rep(1, k - 1)

cocluster_mc(
  n_samples = 1000,
  n_items = n_items,
  target = curve,
  concentration = concentration
)
