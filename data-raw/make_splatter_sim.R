## code to prepare `make_splatter_sim` dataset goes here


## convenience functions

`%>%` <- dplyr::`%>%`
addNAs <- function(group, out.props) rep(NA, nrow(group) - sum(out.props))
addLabels0 <- function(group, out.props, out.sample.names) c(unlist(mapply(rep, times = out.props, x = out.sample.names)), addNAs(group, out.props))

## configs
## the offset is used to generate more cells than required,
## so that we can avoid edge cases
## TODO: setting the seed seems to be problematic
## TODO: need to find a way to make this reproducible

numCellTypes <- 4
numCells <- 500
numGenes <- 2000
numBulkSamples <- 100
numSCsamples <- 10
offset <- 5
total.samples <- numBulkSamples + numSCsamples + offset

seed = 912345
set.seed(seed)

props <- sapply(rep(numCellTypes, total.samples), sample.simplex)
total.props <- round(rowSums(props)/total.samples, 2)
total.props <- total.props/sum(total.props)
sim2 <- splatter::splatSimulate(group.prob = total.props, method = "group", verbose = F,
                                nGenes = numGenes, batchCells = numCells*total.samples, seed = seed)
sim2 <- scater::normalize(sim2)
scater::plotPCA(sim2, colour_by = "Group")

colData.mat <- SingleCellExperiment::colData(sim2)
counts.mat <- SingleCellExperiment::counts(sim2)

## we attach a label signifying to which sample (bulk or single cell) is a particular cell assigned
## we use the colData(sim) DataFrame to keep track of the label

out.props <- round(props[,-((total.samples - offset + 1):total.samples)]*numCells)
out.sample.names <- c(paste0("Bulk", seq(numBulkSamples)), paste0("SC", seq(numSCsamples)))
addLabels <- function(...) addLabels0(out.sample.names = out.sample.names, ...)

## TODO: I can even randomize the labels once they are generated so that we don't have structural problems
groupByCellType <- list()
for (i in seq(numCellTypes)){
  groupByCellType[[i]] <- subset(colData.mat, Group == paste0("Group",i))
  labels <- addLabels(group = groupByCellType[[i]], out.props = out.props[i,])
  groupByCellType[[i]]["label"] <- labels
}
fullSet <- do.call(rbind, groupByCellType)

## for each bulk sample, cells are collapsed to a single sample
splat_bulkExpression <- as.data.frame(matrix(rep(0, numGenes*numBulkSamples), nrow = numGenes, ncol = numBulkSamples))
colnames(splat_bulkExpression) <- paste0("Bulk", seq(numBulkSamples))
for (sample in colnames(splat_bulkExpression)){
  cells <- subset(fullSet, label == sample)$Cell
  gene.exp <- counts.mat[,cells]
  splat_bulkExpression[,sample] <- rowMeans(gene.exp)
}

## for single samples, we collapse the celltypes across samples but retain cellular identity
splat_sigMat <- as.data.frame(matrix(rep(0, numGenes*numCellTypes), nrow = numGenes, ncol = numCellTypes))
colnames(splat_sigMat) <- paste0("Group", seq(numCellTypes))
for (group in colnames(splat_sigMat)){
  cells <- subset(fullSet, Group == group & grepl("SC", label))$Cell
  gene.exp <- counts.mat[,cells]
  splat_sigMat[, group] <- rowMeans(gene.exp)
}

# TODO: I have some testing code here, that needs to be removed.

t(props[,1:5])

out <- baycon(bulkExpression = splat_bulkExpression, sigMat = splat_sigMat, useHyperPrior = T, refresh = 0, iter = 3000)
out$stan$mean

# let's also look at how linear regression does

pEstimate = list()
for (i in 1:5)
  pEstimate[[i]] <- coef(summary(lm(splat_bulkExpression[,i] ~ ., data = splat_sigMat)))[2:5,1]
lmEsts <- do.call(rbind, pEstimate)
lmEsts

# usethis::use_data(splat_bulkExpression, splat_sigMat, overwrite = T)
