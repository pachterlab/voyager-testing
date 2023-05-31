#!/usr/bin/env python
# coding: utf-8

# In[1]:

import voyagerpy as vp
import numpy as np
import pandas as pd
import pathlib
import scanpy as sc
from utils import Checkpoint

checkpoint = Checkpoint()


# ## Quality control (QC)

# In[2]:

# Define the root directory and create it if necessary
root_dir = checkpoint.root_dir / "outs"

adata = vp.read_10x_counts(root_dir, raw=False, datatype="mtx")
is_mt = adata.var["symbol"].str.contains("^mt-", case=False).values
vp.utils.add_per_cell_qcmetrics(adata, subsets={"mito": is_mt})

qc_metrics = [
    "sum",
    "detected",
    "subsets_mito_sum",
    "subsets_mito_detected",
    "subsets_mito_percent",
]
# CHECKPOINT: qc_metrics.csv
checkpoint.add("qc_metrics.csv", adata.obs[qc_metrics])


# In[3]:


qc_features = ["sum", "detected", "subsets_mito_percent"]

cells_to_keep = adata.obs["subsets_mito_percent"] < 20

_, genes_to_keep = np.where(adata[cells_to_keep, :].X.sum(axis=0) > 0)
adata = adata[cells_to_keep, genes_to_keep].copy()

# CHECKPOINT: counts.mtx
checkpoint.add("counts.mtx", adata.X)
# CHECKPOINT: barcodes.txt
checkpoint.add("barcodes.txt", adata.obs_names)
# CHECKPOINT: genes.txt
checkpoint.add("genes.txt", adata.var_names)


# ## Basic non-spatial analysis

# In[4]:


adata.X = adata.X.astype("float64")
adata.layers["counts"] = adata.X.copy()

vp.utils.log_norm_counts(adata, inplace=True)
adata.layers["logcounts"] = adata.X.copy()
# CHECKPOINT: logcounts.mtx
checkpoint.add("logcounts.mtx", adata.layers["logcounts"])


# We use the highly variable genes for PCA.

# In[ ]:


gene_var = vp.utils.model_gene_var(
    adata.layers["logcounts"], gene_names=adata.var_names
)

# CHECKPOINT (sync): gene_var.csv
gene_var = checkpoint.sync("gene_var.csv", gene_var, index_col=0)


# In[ ]:


hvgs = vp.utils.get_top_hvgs(gene_var)
adata.var["highly_variable"] = False
adata.var.loc[hvgs, "highly_variable"] = True
# CHECKPOINT: hvgs.txt
checkpoint.add("hvgs.txt", hvgs)


# In[ ]:


adata.X = vp.utils.scale(adata.X, center=True)
sc.tl.pca(adata, use_highly_variable=True, n_comps=30, random_state=1337)
adata.X = adata.layers["logcounts"].copy()
# CHECKPOINT: pca.mtx
checkpoint.add("pca_embedding.mtx", adata.obsm["X_pca"])
checkpoint.add("pca_vec.mtx", adata[:, hvgs].varm["PCs"])

# In[ ]:


sc.pp.neighbors(
    adata,
    n_pcs=10,
    use_rep="X_pca",
    knn=True,
    n_neighbors=11,
)
sc.tl.leiden(
    adata,
    random_state=29,
    resolution=0.2,
    key_added="cluster",
)

# CHECKPOINT (sync): cluster.csv
clusters = (
    (adata.obs["cluster"].astype(int) + 1).astype(str).astype("category").tolist()
)
adata.obs["clusters"] = checkpoint.sync("cluster.txt", clusters)
adata.obs["clusters"] = (
    (adata.obs["clusters"].astype(int) - 1).astype(str).astype("category")
)


# In[ ]:

markers = vp.utils.find_markers(adata)

# CHECKPOINT (sync): markers_*.csv
for key in markers.keys():
    name, num = key.split("_")
    marker_file: str = f"markers_{int(num)+1}.csv"
    columns = ["p.value", "FDR", "summary.AUC"]
    # columns = ["p_vals", "FDR", "summary_es"]
    df = checkpoint.sync(
        marker_file,
        markers[key],  # .rename(columns=dict(zip(py_columns, r_columns))),
        index_col=0,
    )[columns]
    # df = df.rename(columns=dict(zip(r_columns, py_columns)))
    markers[key] = df.sort_values(by="p.value")


# In[ ]:


marker_genes = [
    marker.sort_values(by="p.value").iloc[0].name
    for _, marker in sorted(markers.items())
]

marker_genes_symbols = adata.var.loc[marker_genes, "symbol"].tolist()


# ## "Spatial" analyses for QC metrics

# In[ ]:


sc.pp.neighbors(
    adata,
    n_neighbors=11,
    n_pcs=10,
    use_rep="X_pca",
    knn=True,
    random_state=29,
    method="gauss",  # one of umap, gauss, rapids
    metric="euclidean",  # many metrics available,
    key_added="knn",
)

dist = adata.obsp["knn_distances"].copy()
dist.data = 1 / dist.data

# row normalize the matrix
dist /= dist.sum(axis=1)
from scipy.sparse import csr_matrix

adata.obsp["knn_weights"] = csr_matrix(dist)

del dist

knn_graph = "knn_weights"

# CHECKPOINT (sync): knn10.mtx
G = checkpoint.sync("knn10.mtx", adata.obsp[knn_graph])
adata.obsp[knn_graph] = G

# adata.obsp["knn_connectivities"] represent the edges, while adata.opsp["knn_weights"] represent the weights
adata.obsp["knn_connectivities"] = (adata.obsp[knn_graph] > 0).astype(int)
vp.spatial.set_default_graph(adata, knn_graph)
_ = vp.spatial.to_spatial_weights(adata, graph_name=knn_graph)


# ### Moran's I

# In[ ]:


vp.spatial.moran(adata, qc_features, graph_name=knn_graph)
morans = adata.uns["spatial"]["moran"][knn_graph].loc[qc_features, ["I"]]

# CHECKPOINT (sync): moran.csv
checkpoint.add("morans.csv", morans.rename(columns={"I": "moran"}), type="csv")


# ### Moran Plot

# In[ ]:


vp.spatial.compute_spatial_lag(adata, qc_features, graph_name=knn_graph, inplace=True)

# CHECKPOINT: moran_plot.csv
lag_features = [f"lagged_{feature}" for feature in qc_features]
df = adata.obs[lag_features].rename(columns=dict(zip(lag_features, qc_features)))
checkpoint.add("moran_plot.csv", df)


# ### Local Moran's I

# Let's compute the locan Moran's metrics for the features.

# In[ ]:


# CHECKPOINT: localmoran.csv
_ = vp.spatial.local_moran(adata, qc_features, graph_name=knn_graph)
n = adata.n_obs
correction = n / (n - 1)
checkpoint.add("localmoran.csv", adata.obsm["local_moran"][qc_features] * correction)


# In[ ]:


# Get the x and y values from PCA
pca_0, pca_1 = adata.obsm["X_pca"][:, :2].T

# Convert the vectors to a gpd.GeoSeries of shapely.Point objects
pca_points = vp.spatial.to_points(pca_0, pca_1)

# Set the points as the geometry for the barcodes
_ = vp.spatial.set_geometry(adata, "pca", pca_points)


# ### Local Spatial heteroscedasticity (LOSH)
# LOSH indicates heterogeneity around each cell in the KNN graph.

# In[ ]:


# CHECKPOINT: losh.csv
_ = vp.spatial.losh(adata, qc_features, graph_name=knn_graph)
checkpoint.add("losh.csv", adata.obsm["losh"][qc_features])


# ## "Spatial" analyses for gene expression

# In[ ]:


top_markers_df = pd.concat(
    [
        pd.DataFrame({"cluster": str(i_cluster), **dict(m[m["FDR"] < 0.05])})
        for i_cluster, (key, m) in enumerate(sorted(markers.items()))
    ]
)

top_markers_df["symbol"] = adata.var.loc[top_markers_df.index, "symbol"]
# top_markers_df


# ### Moran's I

# In[ ]:


# Compute Moran's I for the expression
vp.spatial.moran(adata, feature=hvgs, dim="var", graph_name=knn_graph)
hvgs_moransI = adata.uns["spatial"]["moran"][knn_graph].loc[hvgs, "I"]
# Store the results in adata.var
adata.var.loc[hvgs, "moran"] = hvgs_moransI
# CHECKPOINT: morans_hvgs.csv
checkpoint.add("morans_hvgs.csv", adata.var.loc[hvgs, ["moran"]])


# In[ ]:


# Check if markers any markers are shared across clusters
any(top_markers_df.index.duplicated())


# In[ ]:


adata.var["moran"] = adata.var["moran"].astype("double")

top_markers_df["log_p_adj"] = -np.log10(top_markers_df["FDR"])
top_markers_df["cluster"] = top_markers_df["cluster"].astype("category")
top_markers_df["moran"] = adata.var["moran"]


# In[ ]:


vp.spatial.moran(adata, marker_genes, graph_name=knn_graph, dim="var", permutations=200)
df = adata.uns["spatial"]["moran"]["knn_weights"].loc[marker_genes, ["I"]]
# CHECKPOINT: moran_mc.csv
# TODO: Do we have to add this?
# checkpoint.add("moran_mc.csv", df)


# In[ ]:


vp.spatial.compute_higher_order_neighbors(adata, order=6)
corgram = vp.spatial.compute_correlogram(adata, marker_genes, order=6)
# CHECKPOINT: correlogram.csv
# checkpoint.add("correlogram.csv", corgram)


# ### Local Moran's I

# In[ ]:


_ = vp.spatial.local_moran(adata, marker_genes, layer="logcounts")
# CHECKPOINT:
n = adata.n_obs
correction = n / (n - 1)
checkpoint.add("localmoran.csv", adata.obsm["local_moran"][marker_genes] * correction)


# ### LOSH

# In[ ]:


_ = vp.spatial.losh(adata, marker_genes, graph_name=knn_graph, layer="logcounts")
# CHECKPOINT: losh.csv
checkpoint.add("losh.csv", adata.obsm["losh"][marker_genes])


# ### Moran Plot

# Here we make Moran plots for the top marker genes. We begin by computing the spatial lag explicitly.

# In[ ]:


_ = vp.spatial.compute_spatial_lag(
    adata,
    marker_genes,
    graph_name=knn_graph,
    inplace=True,
    layer="logcounts",
)

lagged_marker_genes = [f"lagged_{feature}" for feature in marker_genes]
df = adata.obs[lagged_marker_genes].rename(
    columns=dict(zip(lagged_marker_genes, marker_genes))
)
# CHECKPOINT: moran_plot.csv
checkpoint.add("moran_plot.csv", df)
