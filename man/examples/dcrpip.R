n_items <- 20
k <- 4
discount <- 0.2

set.seed(42)
samples <- rcrpip(5, n_items, k, discount)
print(samples)

probs <- dcrpip(samples, n_items, k, discount)
print(probs)

log_probs <- dcrpip(samples, n_items, k, discount, log = TRUE)
print(log_probs)

all.equal(exp(log_probs), probs)

single_partition <- c(2, 4, 6, 8)
prob_single <- dcrpip(single_partition, n_items, k, discount)
print(prob_single)

invalid_partition <- c(2, 4, 6, 7)
prob_invalid <- dcrpip(invalid_partition, n_items, k, discount)
print(prob_invalid)
