#!/usr/bin/env Rscript

suppressMessages({
    library(Voyager)
    library(SpatialFeatureExperiment)
    library(SpatialExperiment)
    library(DropletUtils)
    library(BiocNeighbors)
    library(scater)
    library(scran)
    library(bluster)
    library(BiocParallel)
    library(scuttle)
    library(stringr)
    library(BiocSingular)
    library(spdep)
    library(patchwork)
    library(dplyr)
    library(reticulate)
})

root_path <- normalizePath(Sys.getenv("RPATH"))
utils_path <- file.path(root_path, "utils.R")
source(utils_path)
init.dirs()

# Introduction
sce <- read10xCounts(file.path(.root.dir, "outs", "filtered_feature_bc_matrix"))
colnames(sce) <- sce$Barcode
# stop("just checking")

# Quality Control (QC)
is_mito <- str_detect(rowData(sce)$Symbol, "^MT-")
sce <- addPerCellQCMetrics(sce, subsets = list(mito = is_mito))
# CHECKPOINT: qc_metrics.csv
df <- colData(sce)[, c(
	"sum", 
	"detected", 
	"subsets_mito_sum", 
	"subsets_mito_detected", 
	"subsets_mito_percent"
)] %>% as.data.frame
checkpoint.csv("qc_metrics.csv", df)

sce <- sce[, sce$subsets_mito_percent < 20]
sce <- sce[rowSums(counts(sce)) > 0,]


# CHECKPOINT: counts.mtx
. <- checkpoint.mtx( "counts.mtx", counts(sce))
# CHECKPOINT: barcodes.txt
checkpoint.txt("barcodes.txt", colnames(sce))
# CHECKPOINT: genes.txt
checkpoint.txt("genes.txt", rownames(sce))

# Basic non-spatial analyses
sce <- logNormCounts(sce)
# CHECKPOINT: logcounts.mtx
. <- checkpoint.mtx("logcounts.mtx", logcounts(sce))

dec <- modelGeneVar(sce, lowess = FALSE)
hvgs <- getTopHVGs(dec, n = 2000)

# CHECKPOINT (sync): gene_var.csv
checkpoint.csv("gene_var.csv", dec, sync=TRUE)
# CHECKPOINT: hvgs.txt
checkpoint.txt("hvgs.txt", hvgs)

set.seed(29)
sce <- runPCA(sce, ncomponents = 30, BSPARAM = IrlbaParam(),
              subset_row = hvgs, scale = TRUE)

# CHECKPOINT: pca.mtx
pca.mat <- reducedDim(sce, "PCA") %>% Matrix::Matrix(sparse=TRUE)
. <- checkpoint.mtx("pca_embedding.mtx", pca.mat)

pca.rot <- attr(reducedDim(sce, "PCA"), "rotation") %>% Matrix::Matrix(sparse=TRUE)
. <- checkpoint.mtx("pca_vec.mtx", pca.rot)

sce$cluster <- clusterRows(
    reducedDim(sce, "PCA")[,1:10],
    BLUSPARAM = SNNGraphParam(
        cluster.fun = "leiden",
        k = 10,
        cluster.args = list(
            resolution=0.5,
            objective_function = "modularity"
)))

# CHECKPOINT: cluster.txt
checkpoint.txt("cluster.txt", sce$cluster, sync=TRUE)

markers <- findMarkers(sce, groups = colData(sce)$cluster,
                       test.type = "wilcox", pval.type = "all", direction = "up")

# CHECKPOINT: markers.csv
for (i in seq_along(markers)) {
    df <- markers[[i]][c("p.value", "FDR", "summary.AUC")]
    checkpoint.csv(paste0("markers_", i, ".csv"), df, sync=TRUE)
}

top_markers <- unlist(lapply(markers, function(x) head(rownames(x), 1)))
top_markers_symbol <- rowData(sce)[top_markers, "Symbol"]

# "Spatial" analyses for QC metrics

foo <- findKNN(reducedDim(sce, "PCA")[,1:10], k=10, BNPARAM=AnnoyParam())
# Split by row
foo_nb <- asplit(foo$index, 1)
dmat <- 1/foo$distance
# Row normalize the weights
dmat <- sweep(dmat, 1, rowSums(dmat), FUN = "/")
glist <- asplit(dmat, 1)
# Sort based on index
ord <- lapply(foo_nb, order)
foo_nb <- lapply(seq_along(foo_nb), function(i) foo_nb[[i]][ord[[i]]])
class(foo_nb) <- "nb"
glist <- lapply(seq_along(glist), function(i) glist[[i]][ord[[i]]])

listw <- list(style = "W",
              neighbours = foo_nb,
              weights = glist)
class(listw) <- "listw"
attr(listw, "region.id") <- colnames(sce)

# CHECKPOINT: knn10.mtx
. <- checkpoint.knn("knn10.mtx", listw, sync=TRUE)

sfe <- toSpatialFeatureExperiment(sce, spatialCoords = reducedDim(sce, "PCA")[,1:2], spatialCoordsNames = NULL)
colGraph(sfe, "knn10") <- listw

## Moran's I

sfe <- colDataMoransI(sfe, c("sum", "detected", "subsets_mito_percent"))

# CHECKPOINT: morans.csv
moran.cf <- colFeatureData(sfe)[c("sum", "detected", "subsets_mito_percent"), ]["moran_sample01"]
checkpoint.csv( "morans.csv", moran.cf %>% S4Vectors::rename(moran_sample01="moran"))

### Moran plot
sfe <- colDataUnivariate(sfe, "moran.plot", c("sum", "detected", "subsets_mito_percent"))

# CHECKPOINT: moran_plot.csv
lr <- localResults(sfe, name="moran.plot", features=c("sum", "detected", "subsets_mito_percent"))
mp.df <- data.frame(sum=lr[[1]]$wx, detected=lr[[2]]$wx, subsets_mito_percent=lr[[3]]$wx, row.names=rownames(lr[[1]]))
checkpoint.csv("moran_plot.csv", mp.df)
lr <- NULL

### Local Moran's I
sfe <- colDataUnivariate(sfe, "localmoran", c("sum", "detected", "subsets_mito_percent"))

# CHECKPOINT: localmoran.csv
lr <- localResults(sfe, name="localmoran", features=c("sum", "detected", "subsets_mito_percent"))
lm.df <- data.frame(sum=lr[[1]]$Ii, detected=lr[[2]]$Ii, subsets_mito_percent=lr[[3]]$Ii, row.names=rownames(lr[[1]]))
checkpoint.csv("localmoran.csv", lm.df)
lr <- NULL

### Local spatial heteroscadasticity (LOSH)
sfe <- colDataUnivariate(sfe, "LOSH", c("sum", "detected", "subsets_mito_percent"))
# CHECKPOINT: losh.csv
lr <- localResults(sfe, name="LOSH", features=c("sum", "detected", "subsets_mito_percent")) %>% 
    lapply(function(X)as.data.frame(X))
losh.df <- data.frame(sum=lr[[1]]$Hi, detected=lr[[2]]$Hi, subsets_mito_percent=lr[[3]]$Hi, row.names=rownames(lr[[1]]))
rownames(losh.df) <- colnames(sfe)
checkpoint.csv("losh.csv", losh.df)
lr <- NULL

# "Spatial" analyses for gene expression

top_markers_df <- lapply(seq_along(markers), function(i) {
    out <- markers[[i]][markers[[i]]$FDR < 0.05, c("FDR", "summary.AUC")]
    if (nrow(out)) out$cluster <- i
    out
})
top_markers_df <- do.call(rbind, top_markers_df)
top_markers_df$symbol <- rowData(sce)[rownames(top_markers_df), "Symbol"]

## Moran's I
sfe <- runMoransI(sfe, features = hvgs, BPPARAM = MulticoreParam(2))
# CHECKPOINT: morans_hvgs.csv
keep.rows <- !is.na(rowData(sfe)$moran_sample01)
morans <- rowData(sfe)[hvgs, ]["moran_sample01"]
checkpoint.csv("morans_hvgs.csv", morans %>% S4Vectors::rename(moran_sample01="moran"))

top_moran <- head(rownames(sfe)[order(rowData(sfe)$moran_sample01, decreasing = TRUE)], 4)
top_moran_symbol <- rowData(sfe)[top_moran, "Symbol"]
# See if markers are unique to clusters
anyDuplicated(rownames(top_markers_df))

top_markers_df$moran <- rowData(sfe)[rownames(top_markers_df), "moran_sample01"]
top_markers_df$log_p_adj <- -log10(top_markers_df$FDR)
top_markers_df$cluster <- factor(top_markers_df$cluster, 
                                 levels = seq_len(length(unique(top_markers_df$cluster))))

sfe <- runUnivariate(sfe, "moran.mc", features = top_markers, nsim = 200)
# CHECKPOINT: moran_mc.csv
# TODO: Do we have to, we'd only compare the I's?

# TODO: do we need to the correlogram?
# sfe <- runUnivariate(sfe, "sp.correlogram", top_markers, order = 6, 
#                      zero.policy = TRUE, BPPARAM = MulticoreParam(2))

# CHECKPOINT: correlogram.csv
# Not for now.

sfe <- runUnivariate(sfe, "localmoran", features = top_markers)
# CHECKPOINT: localmoran.csv
local_moran_markers <- localResults(sfe, name="localmoran", features=top_markers)

LR <- localResults(sfe, name="localmoran", features=top_markers)
column_names <- paste(top_markers, "Ii", sep=".")
LR <- do.call(cbind, LR)[column_names]
data.table::setnames(LR, old=column_names, new=top_markers)
checkpoint.csv("localmoran.csv", LR)


new_colname <- paste0("cluster", seq_along(top_markers), "_", 
                      top_markers_symbol, "_localmoran")
for (i in seq_along(top_markers)) {
    g <- top_markers[i]
    colData(sfe)[[new_colname[i]]] <- 
        localResult(sfe, "localmoran", g)[,"Ii"]
}

sfe <- runUnivariate(sfe, "LOSH", top_markers)
LR <- localResults(sfe, name="LOSH", features=top_markers) %>% lapply(as.data.frame)
column_names <- paste(top_markers, "Hi", sep=".")
LR <- do.call(cbind, LR)[column_names]
data.table::setnames(LR, old=column_names, new=top_markers)
rownames(LR) <- colnames(sfe)
checkpoint.csv("losh.csv", LR)

new_colname2 <- paste0("cluster", seq_along(top_markers), "_", 
                      top_markers_symbol, "_losh")
                    
for (i in seq_along(top_markers)) {
    g <- top_markers[i]
    colData(sfe)[[new_colname2[i]]] <- 
        localResult(sfe, "LOSH", g)[,"Hi"]
}

sfe <- runUnivariate(sfe, "moran.plot", features = top_markers, colGraphName = "knn10")
# CHECKPOINT: moran_plot.csv
lr <- localResults(sfe, name="moran.plot", features=top_markers)
colnames <- paste(top_markers, "wx", sep=".")
df <- do.call(cbind, lr)[colnames]
data.table::setnames(df, old=colnames, new=top_markers)
checkpoint.csv("moran_plot.csv", df)
