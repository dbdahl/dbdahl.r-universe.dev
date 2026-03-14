# Enumerate all ways to partition 4 items into 2 clusters
partitions <- enumerate_partitions(4, 2)
nrow(partitions)  # S(4,2) = 7
ncol(partitions)  # 4 items

# First partition: items 1,2,3 in cluster 1, item 4 in cluster 2
partitions[1, ]

# Enumerate all ways to partition 3 items into 3 clusters (each item in its own cluster)
partitions <- enumerate_partitions(3, 3)
nrow(partitions)  # S(3,3) = 1
partitions[1, ]   # [1, 2, 3]

# Enumerate all partitions of 4 items (any number of clusters)
partitions <- enumerate_partitions(4)
nrow(partitions)  # B(4) = 15

# Number of set partitions grows rapidly
nrow(enumerate_partitions(5))     # B(5) = 52
nrow(enumerate_partitions(6))     # B(6) = 203
nrow(enumerate_partitions(7))     # B(7) = 877

# Use with n_partitions to verify counts
n_partitions(5)  # B(5) = 52
n_partitions(5, 2)  # S(5,2) = 15
nrow(enumerate_partitions(5, 2))  # Should also be 15
