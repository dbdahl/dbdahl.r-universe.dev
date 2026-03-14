subset <- c("Maine","Georgia","California","Minnesota","Montana")
d <- as.matrix(dist(scale(USArrests[subset,])))
temperature <- 1.0
similarity <- exp( -temperature * d )
DDCRPPartition(similarity=similarity, concentration=0.4)
