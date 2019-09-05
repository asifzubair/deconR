## code to prepare `make_splatter_sim` dataset goes here

## convenience functions
addNAs <- function(group, out.props) rep(NA, nrow(group) - sum(out.props))
addLabels0 <- function(group, out.props, out.sample.names) c(unlist(mapply(rep, times = out.props, x = out.sample.names)),
                                                             addNAs(group, out.props))

## configs
## the offset is used to generate more cells than required,
## so that we can avoid edge cases
## TODO: setting the seed seems to be problematic
## TODO: need to find a way to make this reproducible
numCellTypes <- 4
numCells <- 100
numGenes <- 2000
numBulkSamples <- 100
numSCsamples <- 10
offset <- 5
total.samples <- numBulkSamples + numSCsamples + offset

seed = 912345
set.seed(seed)

splat_props <- sapply(rep(numCellTypes, total.samples), deconR:::sample.simplex)
total.props <- round(rowMeans(splat_props), 2)
total.props <- total.props/sum(total.props)

## Note that if the simplex RV generation works correctly,
## I should have equal proportions in total.props,
## especially when I generate a large number of samples.
splat_sim <- splatter::splatSimulate(group.prob = total.props, method = "group", verbose = F,
                                nGenes = numGenes, batchCells = numCells*total.samples, seed = seed,
                                de.prob = c(0.3, 0.1, 0.2, 0.1),
                                de.downProb = c(0.1, 0.4, 0.6, 0.5),
                                de.facLoc = c(0.6, 0.1, 0.1, 0.2),
                                de.facScale = c(0.1, 0.4, 0.2, 0.5))
splat_sim <- scater::normalize(splat_sim)
## scater::plotPCA(splat_sim, colour_by = "Group")

colData.mat <- SingleCellExperiment::colData(splat_sim)
counts.mat <- SingleCellExperiment::counts(splat_sim)

## we attach a label signifying to which sample (bulk or single cell) is a particular cell assigned
## we use the colData(splat_sim) DataFrame to keep track of the label

remove.offset.samples <- -((total.samples - offset + 1):total.samples)
out.props <- round(splat_props[, remove.offset.samples]*numCells)
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
rownames(splat_bulkExpression) <- rownames(counts.mat)
for (sample in colnames(splat_bulkExpression)){
  cells <- as.character(subset(fullSet, label == sample)$Cell)
  gene.exp <- counts.mat[,cells]
  splat_bulkExpression[,sample] <- rowMeans(gene.exp)
}

## for single cell samples, we collapse the celltypes across samples but retain cellular identity
splat_sigMat <- as.data.frame(matrix(rep(0, numGenes*numCellTypes), nrow = numGenes, ncol = numCellTypes))
colnames(splat_sigMat) <- paste0("Group", seq(numCellTypes))
rownames(splat_sigMat) <- rownames(counts.mat)
for (group in colnames(splat_sigMat)){
  cells <- as.character(subset(fullSet, Group == group & grepl("SC", label))$Cell)
  gene.exp <- counts.mat[,cells]
  splat_sigMat[, group] <- rowMeans(gene.exp)
}

## save it
remove.sc.offset.samples <- -((total.samples - numSCsamples - offset + 1):total.samples)
splat_props <- splat_props[ , remove.sc.offset.samples]
usethis::use_data(splat_props, splat_sigMat, splat_bulkExpression, overwrite = T)
