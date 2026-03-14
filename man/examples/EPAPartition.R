subset <- c("Maine","Georgia","California","Minnesota","Montana")
d <- as.matrix(dist(scale(USArrests[subset,])))
temperature <- 1.0
similarity <- exp( -temperature * d )
EPAPartition(similarity=similarity, permutation=c(1,5,4,2,3), concentration=1.5, discount=0.1)
