# Partition n=20 items into k=4 clusters with target weights 1:2:3:4
n_items <- 20
weights <- c(1, 2, 3, 4) # automatically sorted and normalized internally

# TiDaL kernel (default): omit tail_shape
samples_tidal <- rlorenzip(500, n_items, weights, concentration = 10)
print(samples_tidal[1:5, ])

# TaDPoLe kernel: supply tail_shape
samples_tadpole <- rlorenzip(
  500, n_items, weights,
  concentration = 10,
  tail_shape = 0.5
)
print(samples_tadpole[1:5, ])

# Using a Lorenz curve object as target (n_clusters required)
lc <- lorenz_linear(c(1, 2, 3, 7))
samples_lc <- rlorenzip(100, n_items, lc, n_clusters = 4, concentration = 5)
print(samples_lc[1:5, ])
