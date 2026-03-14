# Integer partition function p(n): total number of integer partitions
n_integer_partitions(0) # 1 (empty partition)
n_integer_partitions(1) # 1
n_integer_partitions(5) # 7 (5, 1:4, 2:3, 1:1:3, 1:2:2, 1:1:1:2, 1:1:1:1:1)
n_integer_partitions(10) # 42

# Number of integer partitions into exactly k parts: p(n, k)
n_integer_partitions(5, n_clusters = 1) # 1 (just 5)
n_integer_partitions(5, n_clusters = 2) # 2 (1:4, 2:3)
n_integer_partitions(5, n_clusters = 3) # 2 (1:1:3, 1:2:2)
n_integer_partitions(5, n_clusters = 4) # 1 (1:1:1:2)
n_integer_partitions(5, n_clusters = 5) # 1 (1:1:1:1:1)

# Verify: sum of p(n, k) over k equals p(n)
sum(sapply(1:5, function(k) n_integer_partitions(5, k))) # 7 = p(5)

# Verify against enumerate_integer_partitions
n_items <- 10
n_clusters <- 3
nrow(enumerate_integer_partitions(n_items, n_clusters)) # number of rows
n_integer_partitions(n_items, n_clusters) # should match

# Use log = TRUE for large values to avoid overflow
n_integer_partitions(100) # 190569292
n_integer_partitions(200, log = TRUE) # log of p(200)
exp(n_integer_partitions(200, log = TRUE)) # approximately 3.97e12

# Compare integer partitions vs set partitions (Bell numbers)
# These are different concepts!
n_integer_partitions(5) # 7 (ways to write 5 as sum of positive integers)
n_partitions(5) # 52 (ways to set-partition 5 distinguishable items)
