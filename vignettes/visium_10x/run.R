#!/usr/bin/env Rscript

suppressMessages({
	library(Voyager)
	library(SpatialExperiment)
	library(SpatialFeatureExperiment)
	library(SingleCellExperiment)
	library(ggplot2)
	library(scater)
	library(scuttle)
	library(scran)
	library(stringr)
	library(patchwork)
	library(bluster)
	library(rjson)
})

root_path <- normalizePath(Sys.getenv("RPATH"))
utils_path <- file.path(root_path, "utils.R")
source(utils_path)
init.dirs()

sfe <- read10xVisiumSFE(samples = .root.dir, type = "sparse", data = "raw")
is_mt <- str_detect(rowData(sfe)$symbol, "^mt-")
sfe <- addPerCellQCMetrics(sfe, subsets = list(mito = is_mt))

# CHECKPOINT: qc_metrics.csv
df <- colData(sfe)[, c(
	"sum", 
	"detected", 
	"subsets_mito_sum", 
	"subsets_mito_detected", 
	"subsets_mito_percent"
)] %>% as.data.frame
checkpoint.csv("qc_metrics.csv", df)

sfe <- SpatialFeatureExperiment::transpose(sfe)

sfe_tissue <- sfe[,sfe$in_tissue]
# CHECKPOINT: tissue_ids.txt
checkpoint.txt("tissue_ids.txt", colnames(sfe_tissue))

# CHECKPOINT: counts.mtx
. <- checkpoint.mtx("counts.mtx", counts(sfe_tissue))

sfe_tissue <- logNormCounts(sfe_tissue)
# CHECKPOINT: logcounts.mtx
. <- checkpoint.mtx("logcounts.mtx", logcounts(sfe_tissue))

dec <- modelGeneVar(sfe_tissue, lowess = FALSE)
# CHECKPOINT (sync): gene_var.csv
checkpoint.csv("gene_var.csv", dec, sync=TRUE)

hvgs <- getTopHVGs(dec, n = 2000)

# CHECKPOINT: hvgs.txt
checkpoint.txt("hvgs.txt", hvgs, sync = FALSE) # we want to compare hvgs given gene_var
set.seed(29)
sfe_tissue <- runPCA(sfe_tissue, ncomponents = 30, subset_row = hvgs,
                     scale = TRUE, BSPARAM=BiocSingular::ExactParam()) # scale as in Seurat

# CHECKPOINT: pca.mtx

pca.mat <- reducedDim(sfe_tissue, "PCA") %>% Matrix::Matrix(sparse=TRUE)
. <- checkpoint.mtx("pca.mtx", pca.mat)

set.seed(29)
colData(sfe_tissue)$cluster <- clusterRows(
	reducedDim(sfe_tissue, "PCA")[,1:3],
	BLUSPARAM = SNNGraphParam(
		cluster.fun = "leiden",
		cluster.args = list(
			resolution_parameter = 0.5,
			objective_function = "modularity"
		)
	)
)

# CHECKPOINT: cluster.txt
checkpoint.txt("cluster.txt", sfe_tissue$cluster, sync=TRUE)

markers <- findMarkers(sfe_tissue, groups = colData(sfe_tissue)$cluster,
                       test.type = "wilcox", pval.type = "all", direction = "up")

# TODO: Write the markers to file(s) for each cluster

genes_use <- vapply(markers, function(x) rownames(x)[1], FUN.VALUE = character(1))
# CHECKPOINT (sync): genes_use.txt
checkpoint.txt("marker_genes_set.txt", genes_use, sync=TRUE)