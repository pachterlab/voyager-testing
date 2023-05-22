#!/usr/bin/env python

import numpy as np
import geopandas as gpd

import scanpy as sc
import voyagerpy as vp

import pathlib

from utils import Checkpoint

checkpoint = Checkpoint()

outs_dir = pathlib.Path(checkpoint.root_dir / "outs")

adata = vp.read_10x_visium(
    outs_dir,
    datatype="mtx",
    raw=True,
    prefix=None,
    symbol_as_index=False,
    dtype=np.float64,
    res="lowres",
)

is_mt = adata.var["symbol"].str.contains("^mt-").values
vp.utils.add_per_cell_qcmetrics(adata, subsets={"mito": is_mt})

df = adata.obs[
    [
        "sum",
        "detected",
        "subsets_mito_sum",
        "subsets_mito_detected",
        "subsets_mito_percent",
    ]
]

# CHECKPOINT: qc_metrics.csv
checkpoint.add("qc_metrics.csv", df)

scale = vp.utils.get_scale(adata)
spot_diam = adata.uns["spatial"]["scale"]["spot_diameter_fullres"]

# Create the polygons to represent the Visium spots
visium_spots = vp.spatial.to_points(
    x="pxl_col_in_fullres",
    y="pxl_row_in_fullres",
    data=adata.obs,
    scale=scale,
    radius=scale * spot_diam / 2,
)

# or, equivalently
spots_x, spots_y = vp.spatial.get_spot_coords(adata)
spots = gpd.GeoSeries.from_xy(spots_x, spots_y, index=adata.obs_names).buffer(
    scale * spot_diam / 2
)

del spots

# Set the geometry to the visium spots and assign the name "spot_poly"
_ = vp.spatial.set_geometry(adata, geom="spot_poly", values=visium_spots)

# set the `in_tissue` as a categorical variable.
adata.obs["in_tissue"] = adata.obs["in_tissue"].astype(bool).astype("category")

adata_tissue = adata[adata.obs["in_tissue"] == True].copy()
vp.spatial.set_geometry(adata_tissue, "spot_poly")

# CHECKPOINT: tissue_ids.txt
checkpoint.add("tissue_ids.txt", adata_tissue.obs_names)

# The original count data
adata_tissue.layers["counts"] = adata_tissue.X.copy()

# CHECKPOINT: counts.mtx
checkpoint.add("counts.mtx", adata_tissue.layers["counts"])

# Log-normalize the adata.X matrix
vp.utils.log_norm_counts(adata_tissue, inplace=True)
adata_tissue.layers["logcounts"] = adata_tissue.X.copy()

# CHECKPOINT: logcounts.mtx
checkpoint.add("logcounts.mtx", adata_tissue.layers["logcounts"])

gene_var = vp.utils.model_gene_var(
    adata_tissue.layers["logcounts"], gene_names=adata_tissue.var_names
)

# CHECKPOINT: gene_var.csv
gene_var = checkpoint.sync("gene_var.csv", gene_var, index_col=0)

# Compute the highly variable genes
hvgs = vp.utils.get_top_hvgs(gene_var)

# CHECKPOINT: hvgs.txt
checkpoint.add("hvgs.txt", hvgs, type="txt")
# sync?

# Set the 'highly_variable' column for the genes
adata_tissue.var["highly_variable"] = False
adata_tissue.var.loc[hvgs, "highly_variable"] = True

# scale first, then perform pca

adata_tissue.X = vp.utils.scale(adata_tissue.X, center=True)
sc.tl.pca(adata_tissue, use_highly_variable=True, n_comps=30, random_state=1337)

# CHECKPOINT: pca.mtx
checkpoint.add("pca.mtx", adata_tissue.obsm["X_pca"])

from leidenalg import ModularityVertexPartition

sc.pp.neighbors(adata_tissue, n_pcs=3, use_rep="X_pca", method="gauss", n_neighbors=80)
sc.tl.leiden(
    adata_tissue,
    random_state=29,
    resolution=None,
    key_added="cluster",
    partition_type=ModularityVertexPartition,
)

# CHECKPOINT: cluster.txt; we add one since R is 1-indexed
clusters = (adata_tissue.obs["cluster"].astype(int) + 1).astype(str)
checkpoint.sync("cluster.txt", clusters.tolist(), type="txt")

markers = vp.utils.get_marker_genes(adata_tissue, False, cluster="cluster")
marker_genes = markers.iloc[0, :].tolist()

# CHECKPOINT: marker_genes.txt;
checkpoint.sync("marker_genes_set.txt", marker_genes, type="txt")
