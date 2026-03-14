# Create a Lorenz curve from an empirical partition
partition <- c(1, 2, 3, 4, 5)

# Piecewise linear Lorenz curve
curve_linear <- lorenz_linear(partition)
print(curve_linear)

# I-spline Lorenz curve (smooth)
curve_ispline <- lorenz_ispline(partition)
print(curve_ispline)

# Evaluate curves at specific points
x <- seq(0, 1, by = 0.25)
lorenz_evaluate(curve_linear, x)
lorenz_evaluate(curve_ispline, x)

# Compute Gini coefficient
gini(curve_linear)
gini(curve_ispline)

# Direct Gini computation from partition (without creating curve)
gini_coefficient(partition)

# Extract Lorenz curve increments for a different number of clusters
lorenz_increments(curve_linear, n_clusters = 3)
lorenz_increments(curve_ispline, n_clusters = 3)

# Use Lorenz curve with rlorenzip
set.seed(42)
samples <- rlorenzip(
  n = 5,
  n_items = 20,
  target = curve_linear,
  n_clusters = 5,
  concentration = 10
)
print(samples)

# Plot a Lorenz curve (if graphics available)
if (interactive()) {
  plot(curve_linear)
  plot(curve_ispline)
}
