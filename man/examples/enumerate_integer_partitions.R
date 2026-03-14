# Enumerate all integer partitions of 10 into 3 parts
partitions <- enumerate_integer_partitions(10, 3)
print(partitions)

# Each row sums to n_items
rowSums(partitions)

# Each row is in non-decreasing order
all(apply(partitions, 1, function(x) all(diff(x) >= 0)))

# Number of integer partitions grows with n_items and n_clusters
nrow(enumerate_integer_partitions(10, 2))
nrow(enumerate_integer_partitions(10, 3))
nrow(enumerate_integer_partitions(10, 4))
nrow(enumerate_integer_partitions(20, 4))

# Use with dlorenzip to compute exact probabilities over all integer partitions
n_items <- 10
n_clusters <- 3
weights <- c(1, 2, 3)
concentration <- 5

all_partitions <- enumerate_integer_partitions(n_items, n_clusters)
probs <- dlorenzip(all_partitions, n_items, weights, concentration = concentration)

# Probabilities sum to 1
sum(probs)

# Find the most likely integer partition
all_partitions[which.max(probs), ]
