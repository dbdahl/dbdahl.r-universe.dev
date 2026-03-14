n_items <- 80
crp_concentration <- 1
crp_discount <- 0.2

target <- lorenz_linear(c(1, 2, 7, 12))
plot(target)
concentration <- 1

set.seed(123)
out <- crp_vs_lorenz(
  n_items,
  target,
  concentration,
  crp_concentration = crp_concentration,
  crp_discount = crp_discount,
  n_samples = 200
)
print(out)
