library(Voyager)
library(SpatialExperiment)

library(SpatialFeatureExperiment)
library(SingleCellExperiment)

library(scater)
#> Loading required package: scuttle
library(scuttle)
library(scran)
library(stringr)
library(patchwork)
library(bluster)
library(rjson)

library(Matrix)


if (!file.exists("visium_ob.tar.gz"))
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_raw_feature_bc_matrix.tar.gz", 
                  destfile = "visium_ob.tar.gz")

if (!dir.exists("outs")) {
    dir.create("outs")
    system("tar -xvf visium_ob.tar.gz -C outs")
    system("tar -xvf visium_ob_spatial.tar.gz -C outs")
}


fromJSON(file = "outs/spatial/scalefactors_json.json")
sfe <- read10xVisiumSFE(samples = "outs", type = "sparse", data = "raw")


is_mt <- str_detect(rowData(sfe)$symbol, "^mt-")
sfe <- addPerCellQCMetrics(sfe, subsets = list(mito = is_mt))

sfe_tissue <- sfe[,sfe$in_tissue]
sfe_tissue <- logNormCounts(sfe_tissue)


lc <- assay(sfe_tissue, 'logcounts')
writeMM(lc, 'checkpoints/checkpoint_1_vr.mtx')


sfe_tissue <- runPCA(sfe_tissue, ncomponents = 30,# subset_row = hvgs,
                     scale = TRUE) # scale as in Seurat

writeMM(reducedDim(sfe_tissue, "PCA"), 'checkpoints/checkpoint_2_vr.mtx')
