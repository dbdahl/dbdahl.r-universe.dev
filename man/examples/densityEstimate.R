example("nealAlgorithm")

density <- function(grid, parameter) {
  dnorm(grid, mean=parameter, sd=sqrt(sigma2))
}

z <- pumpkin:::densityEstimate(partitionSamples, parametersSamples, distr, sampleParameter, density, seq(-3,3,length=1000))
plot(z,type="l", xlab="Response", ylab="Density")
for ( i in seq_along(data) ) {
  lines(c(data[i],data[i]),c(0,0.05))
}
