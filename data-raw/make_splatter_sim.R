## code to prepare `make_splatter_sim` dataset goes here

numCellTypes <- 4
numCells <- 100
numGenes <- 2001
numBulkSamples <- 100
numSCsamples <- 6
total.samples <- numBulkSamples + numSCsamples
seed = 912345

set.seed(seed)
props <- sapply(rep(numCellTypes, total.samples), sample.simplex)
total.props <- round(rowSums(props)/total.samples, 2)
total.props <- total.props/sum(total.props)
sim2 <- splatter::splatSimulate(group.prob = total.props, method = "group", verbose = F,
                                nGenes = numGenes, batchCells = numCells*total.samples, seed = seed)
sim2 <- scater::normalize(sim2)
scater::plotPCA(sim2, colour_by = "Group")



#usethis::use_data()
