# Bell numbers: total number of set partitions of n items
n_partitions(0)
n_partitions(1)
n_partitions(5)
n_partitions(10)

# Stirling numbers of the second kind: S(n, k)
# Number of ways to set-partition n items into exactly k non-empty groups
n_partitions(5, n_clusters = 1) # 1 (all items in one group)
n_partitions(5, n_clusters = 2) # 15
n_partitions(5, n_clusters = 3) # 25
n_partitions(5, n_clusters = 4) # 10
n_partitions(5, n_clusters = 5) # 1 (each item in its own group)

# Verify: sum of Stirling numbers equals Bell number
sum(sapply(1:5, function(k) n_partitions(5, k))) # 52 = B(5)

# Multinomial coefficient for set partitions
# Number of ways to set-partition n items into groups of specific sizes
n_partitions(7, n_clusters = 3, integer_partition = c(2, 2, 3))
n_partitions(6, n_clusters = 3, integer_partition = c(2, 2, 2))
n_partitions(6, n_clusters = 2, integer_partition = c(3, 3))

# Verify: sum over all integer partitions equals Stirling number
n_items <- 10
n_clusters <- 3
all_integer_partitions <- enumerate_integer_partitions(n_items, n_clusters)
total <- sum(apply(all_integer_partitions, 1, function(p) {
  n_partitions(n_items, n_clusters, p)
}))
total # Should equal S(10, 3)
n_partitions(n_items, n_clusters) # S(10, 3) = 9330

# Use log = TRUE for large values to avoid overflow
n_partitions(20) # 51724158235372
n_partitions(25, log = TRUE) # log of B(25)
exp(n_partitions(25, log = TRUE)) # approximately 4.64e18

# Verify relationship with dlorenzip: probabilities sum to 1
n_items <- 10
n_clusters <- 3
weights <- c(1, 2, 3)
concentration <- 5

all_integer_partitions <- enumerate_integer_partitions(n_items, n_clusters)
probs <- dlorenzip(all_integer_partitions, n_items, weights, concentration = concentration)
sum(probs) # Should be close to 1

# Expected number of set partitions for a random integer partition
# drawn from the Lorenz IP distribution
expected_n_set_partitions <- sum(sapply(
  seq_len(nrow(all_integer_partitions)),
  function(i) {
    probs[i] * n_partitions(n_items, n_clusters, all_integer_partitions[i, ])
  }
))
expected_n_set_partitions

# Compare to the average under a uniform distribution over integer partitions
n_partitions(n_items, n_clusters) / nrow(all_integer_partitions)
