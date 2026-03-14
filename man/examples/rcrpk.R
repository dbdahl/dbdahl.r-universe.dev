# Standard CRP (Dirichlet process) with concentration = 1
set.seed(42)
k_samples <- rcrpk(1000, n_items = 100, concentration = 1, discount = 0)
hist(k_samples, main = "Number of clusters (CRP, alpha=1)")
mean(k_samples)  # Expected ~ alpha * log(n)

# Two-parameter CRP with discount
k_samples_py <- rcrpk(1000, n_items = 100, concentration = 1, discount = 0.5)
hist(k_samples_py, main = "Number of clusters (Pitman-Yor)")
mean(k_samples_py)  # Generally more clusters than standard CRP
