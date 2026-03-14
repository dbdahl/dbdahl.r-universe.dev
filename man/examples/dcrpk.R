# PMF for standard CRP
k_values <- 1:20
probs <- dcrpk(k_values, n_items = 50, concentration = 2, discount = 0)
barplot(probs, names.arg = k_values, xlab = "Number of clusters", ylab = "Probability")

# Verify probabilities sum to 1
sum(dcrpk(1:50, n_items = 50, concentration = 2, discount = 0))

# Log-probabilities for numerical stability with large n_items
log_probs <- dcrpk(1:100, n_items = 1000, concentration = 5, discount = 0.25, log = TRUE)

# Compare CRP vs Pitman-Yor
k <- 1:30
p_crp <- dcrpk(k, n_items = 100, concentration = 1, discount = 0)
p_py <- dcrpk(k, n_items = 100, concentration = 1, discount = 0.5)
plot(k, p_crp, type = "h", col = "blue", lwd = 2,
     xlab = "Number of clusters", ylab = "Probability")
lines(k + 0.2, p_py, type = "h", col = "red", lwd = 2)
legend("topright", c("CRP (d=0)", "Pitman-Yor (d=0.5)"), col = c("blue", "red"), lwd = 2)
