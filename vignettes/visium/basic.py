#!/usr/bin/env python
# coding: utf-8



import networkx as nx
import numpy as np
import geopandas as gpd
import pandas as pd

import scanpy as sc
import voyagerpy as vp
import os

from anndata import AnnData
import pathlib
import json
import requests

outs_dir = pathlib.Path('outs')
outs_dir.mkdir(parents=True, exist_ok=True)
root_dir = (outs_dir / '..').resolve()


checkpath = root_dir / 'checkpoints'
checkpath.mkdir(parents=True, exist_ok=True)

# Download the gene count matrix
tar_path_ob = root_dir / 'visium_ob.tar.gz'
url_reads = "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_raw_feature_bc_matrix.tar.gz"
if not tar_path_ob.exists():
    res = requests.get(url_reads)
    with tar_path_ob.open('wb') as f:
        f.write(res.content)

# Download the spatial information
tar_path_sp =  root_dir / 'visium_ob_spatial.tar.gz'
url_spatial = "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz"
if not tar_path_sp.exists():
    res = requests.get(url_spatial)
    with tar_path_sp.open('wb') as f:
        f.write(res.content)

# Decompress the downloaded files, TODO only open if neeeded
os.system(f'tar -xvf {tar_path_ob} -C {outs_dir}')
os.system(f'tar -xvf {tar_path_sp} -C {outs_dir}')



adata = vp.read_10x_visium(
    outs_dir,
    datatype = 'mtx',
    raw = True,
    prefix = None,
    symbol_as_index=False
)



vp.plt.imshow(adata)


if adata.uns['config']['var_names'] == 'symbol':
    is_mt = adata.var.index.to_series().str.contains('^mt-').values
else:
    is_mt = adata.var['symbol'].str.contains('^mt-').values

adata = vp.utl.calculate_metrics(adata)
vp.utils.add_per_cell_qcmetrics(adata, subsets={'mito': is_mt})

adata = vp.spatial.get_geom(adata)

#_ = vp.plt.plot_spatial_feature(
##    adata,
#    ['sum', 'detected', 'subsets_mito_percent'],
#    tissue=False,
#    ncol=2
#)


adata.obs['in_tissue'] = adata.obs['in_tissue'].astype('category')

adata_tissue = adata[adata.obs["in_tissue"]==1]
adata_tissue.obs = gpd.GeoDataFrame(adata_tissue.obs, geometry="spot_poly")



adata_tissue.X = adata_tissue.X.astype('float64')
adata_tissue.layers['counts'] = adata_tissue.X.copy()

target_sum = adata_tissue.X.sum(axis=1).mean()
sc.pp.normalize_total(adata_tissue, target_sum=target_sum)
sc.pp.log1p(adata_tissue, base=2)
adata_tissue.layers['logcounts'] = adata_tissue.X.copy()



from scipy import sparse, io


checkpath = root_dir / 'checkpoints'

io.mmwrite(checkpath / 'checkpoint_1_vp.mtx', adata_tissue.X)

#adata_tissue.X = vp.utils.scale(adata_tissue.X, center=True)
sc.tl.pca(adata_tissue, use_highly_variable=False, n_comps=30, random_state=1337)

io.mmwrite(checkpath / 'checkpoint_2_vp.mtx', adata_tissue.obsm['X_pca'])

