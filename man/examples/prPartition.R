concentration <- 1.0
discount <- 0.1
nSamples <- 3

distr <- CRPPartition(nItems=5, concentration=concentration, discount=discount)
x <- samplePartition(distr, nSamples, nCores=1)
prPartition(distr, x)

anchor <- c(1,1,1,2,2)
permutation <- c(1,5,4,2,3)
n_items <- length(permutation)

distr <- LocationScalePartition(anchor=anchor, shrinkage=1/0.1,
                                concentration=concentration, permutation=permutation)
x <- samplePartition(distr, nSamples, nCores=1)
prPartition(distr, x)

distr <- CenteredPartition(anchor=anchor, shrinkage=1/0.1,
             CRPPartition(nItems=n_items, concentration=concentration, discount=discount))
x <- samplePartition(distr, nSamples, nCores=1)
prPartition(distr, x)

distr <- ShrinkagePartition(anchor=anchor, shrinkage=c(0,0,0,0.3,0.3),
             permutation=permutation, grit=0.2,
             CRPPartition(nItems=n_items, concentration=concentration, discount=discount))
x <- samplePartition(distr, nSamples, nCores=1)
prPartition(distr, x)

subset <- c("Maine","Georgia","California","Minnesota","Montana")
d <- as.matrix(dist(scale(USArrests[subset,])))
temperature <- 1.0
similarity <- exp( - temperature * d )
distr <- EPAPartition(similarity=similarity, permutation=permutation,
                      concentration=concentration, discount=discount)
x <- samplePartition(distr, nSamples, nCores=1)
prPartition(distr, x)

distr <- DDCRPPartition(similarity=similarity, concentration=1.5)
samplePartition(distr, nSamples, nCores=1)
