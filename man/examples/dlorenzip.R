n_items <- 20
weights <- c(1, 2, 3, 4)
concentration <- 1

# Sample some integer partitions (target is a numeric vector of weights)
set.seed(42)
samples <- rlorenzip(5, n_items, weights, concentration = concentration)
print(samples)

# Compute probabilities for the sampled integer partitions
probs <- dlorenzip(samples, n_items, weights, concentration = concentration)
print(probs)

# Compute log probabilities
log_probs <- dlorenzip(samples, n_items, weights, concentration = concentration, log = TRUE)
print(log_probs)

# Verify consistency: exp(log_probs) should equal probs
all.equal(exp(log_probs), probs)

# Compute probability for a single integer partition (as a vector)
single_partition <- c(2, 4, 6, 8)
prob_single <- dlorenzip(single_partition, n_items, weights, concentration = concentration)
print(prob_single)

# Using the TaDPoLe variant with tail-shape parameters
log_skew <- 0
tail_shape <- 1.0
probs_tadpole <- dlorenzip(
  samples,
  n_items,
  weights,
  concentration = concentration,
  log_skew = log_skew,
  tail_shape = tail_shape
)
print(probs_tadpole)

# Invalid integer partitions return probability 0 (or -Inf for log)
invalid_partition <- c(2, 4, 6, 7) # sums to 19, not 20
prob_invalid <- dlorenzip(invalid_partition, n_items, weights, concentration = concentration)
print(prob_invalid) # 0

log_prob_invalid <- dlorenzip(
  invalid_partition,
  n_items,
  weights,
  concentration = concentration,
  log = TRUE
)
print(log_prob_invalid) # -Inf

# Using target as a Lorenz curve object
# Use a more unequal integer partition for higher Gini
lorenz_curve <- lorenz_linear(c(1, 2, 3, 7))
probs_from_lorenz <- dlorenzip(
  samples,
  n_items,
  lorenz_curve,
  n_clusters = 4,
  concentration = concentration
)
print(probs_from_lorenz)
