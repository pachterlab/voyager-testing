library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(Matrix)

# Download filtered data from 10x website
if (!file.exists("visium_ob.tar.gz"))
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.tar.gz",
                  destfile = "visium_ob.tar.gz")
if (!file.exists("visium_ob_spatial.tar.gz"))
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz", 
                  destfile = "visium_ob_spatial.tar.gz")
if (!dir.exists("outs")) {
    dir.create("outs")
    system("tar -xvf visium_ob.tar.gz -C outs")
    system("tar -xvf visium_ob_spatial.tar.gz -C outs")
}

sfe <- read10xVisiumSFE(samples = "outs", type = "sparse", data = "filtered")
# Use a small subset of data for testing
sfe_sub <- crop(sfe, xmin = 3000, xmax = 4000, ymin = 4000, ymax = 5000)
sfe_sub <- sfe_sub[rowSums(counts(sfe_sub)) > 10,]
# Randomly select 50 genes
set.seed(29)
genes_use <- sample(rownames(sfe_sub), 50)
sfe_sub <- sfe_sub[genes_use,]
sfe_sub <- logNormCounts(sfe_sub)

# Save the log normalized matrix
dir.create("testdata")
writeMM(logcounts(sfe_sub), file = "testdata/logcounts.mtx")
writeMM(counts(sfe_sub), file = "testdata/counts.mtx")
writeLines(colnames(sfe_sub), con = file("testdata/barcodes.tsv"))
writeLines(rownames(sfe_sub), con = file("testdata/genes.tsv"))
write.csv(spatialCoords(sfe_sub), file = "testdata/coords.csv", quote = FALSE, 
          row.names = FALSE)
st_write(spotPoly(sfe_sub), dsn = "testdata/spotPoly.geojson")
