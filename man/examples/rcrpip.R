# Sample integer partitions from the CRP (Ewens case, discount = 0)
set.seed(42)
samples <- rcrpip(n = 5, n_items = 20, n_clusters = 4, discount = 0)
print(samples)

# With discount > 0 (Pitman-Yor case)
samples_py <- rcrpip(n = 5, n_items = 20, n_clusters = 4, discount = 0.25)
print(samples_py)
