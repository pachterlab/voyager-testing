import numpy as np
import geopandas as gpd
import pandas as pd

import scanpy as sc
import voyagerpy as vp
from scipy import sparse, io
import os

from anndata import AnnData
from pathlib import Path
from typing import List, Optional, Sequence, Union

#%%
import requests

import functools
from collections import defaultdict

checkpoint_counters = defaultdict(int)


def checkpoint_fun(fun):
    @functools.wraps(fun)
    def inner(
        *args,
        filename: Union[str, Path],
        dir: Optional[Path] = None,
        increment: bool = True,
        include_counter: bool = True,
        **kwargs,
    ):
        global checkpoint_counters
        filename = Path(filename)

        if dir is not None:
            filename = dir / filename

        filename = filename.resolve()
        parent_dir = filename.parent
        parent_dir.mkdir(parents=True, exist_ok=True)

        i_checkpoint = checkpoint_counters[parent_dir]

        if include_counter:
            filename = filename.with_stem(f"{filename.stem}_{i_checkpoint}")
        print("Saving checkpoint", filename)
        ret = fun(*args, filename=filename, **kwargs)

        if include_counter and increment:
            checkpoint_counters[parent_dir] += 1
        return ret

    return inner


@checkpoint_fun
def mtx_checkpoint(df, *, filename: Path, **kwargs):
    """write df to mtx file. All elements must be numbers"""
    io.mmwrite(filename, df)


@checkpoint_fun
def txt_checkpoint(
    items: List[Union[str, Sequence[str]]], filename: Path, sep: str = ",", **kwargs
):
    """write items into text file. If an element in items is an iterable of str,
    they are concatenated with sep.

    items: List[Union[str, Sequence[str]]] - the items to write
    filename: Path - the path to the output file
    sep: str - the separator to use if an item in items is a sequence
    """
    with filename.open("wt") as f:
        rows = map(sep.join, items)
        f.writelines([row + "\n" for row in rows])


def main(cwd: Path, checkpoint_dir: Path, data_dir: Path, sync_with_other: bool = True):
    download_data(cwd, data_dir)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    sync_dir = cwd / "sync_data"

    adata = vp.read_10x_visium(
        data_dir,
        datatype="mtx",
        raw=True,
        prefix=None,
        symbol_as_index=False,
        dtype=np.float64,
    )

    # rotate image and coordinates. Should not affect outcome
    vp.spatial.rotate_img90(adata, k=-1, apply=True)
    vp.spatial.mirror_img(adata, axis=1, apply=True)

    # Quality Control (QC)
    if adata.uns["config"]["var_names"] == "symbol":
        is_mt = adata.var.index.to_series().str.contains("^mt-").values
    else:
        is_mt = adata.var["symbol"].str.contains("^mt-").values

    adata = vp.utils.calculate_metrics(adata)
    vp.utils.add_per_cell_qcmetrics(adata, subsets={"mito": is_mt})

    # Is this necessary? Probaby for the spatial analyses
    adata = vp.spatial.get_geom(adata)

    # TODO: add checkpoint before this plot
    checkpoint_features = [
        "sum",
        "detected",
        "subsets_mito_sum",
        "subsets_mito_detected",
        "subsets_mito_percent",
    ]
    mtx_checkpoint(
        adata.obs[checkpoint_features],
        filename=checkpoint_dir / "chkpt.mtx",
    )

    # PLOT: plot_spatial_feature(adata, ['sum', 'detected', 'subsets_mito_percent'])

    adata.obs["in_tissue"] = adata.obs["in_tissue"].astype("category")
    # PLOT: plot_barcode_data as violinplots, same features as above. Don't need a checkpoint
    # PLOT: plot_barcode_data as scatter (subsets_mito_percent and sum). Don't need a checkpoint

    adata_tissue = adata[adata.obs["in_tissue"] == 1]
    del adata

    txt_checkpoint(
        adata_tissue.obs.index.tolist(), filename=checkpoint_dir / "tissue_spots.txt"
    )

    # Do we need this?
    adata_tissue = vp.spatial.get_geom(adata_tissue)

    adata_tissue.layers["counts"] = adata_tissue.X.copy()
    vp.utils.log_norm_counts(adata_tissue, inplace=True)
    adata_tissue.layers["logcounts"] = adata_tissue.X.copy()

    # TODO: Add checkpoint to verify subset: counts layer vs counts assay
    mtx_checkpoint(
        adata_tissue.layers["counts"], filename=checkpoint_dir / "counts.mtx"
    )
    # TODO: Add checkpoint to verify log_norm_counts: logcounts layer vs logcounts assay
    mtx_checkpoint(
        adata_tissue.layers["logcounts"], filename=checkpoint_dir / "logcounts.mtx"
    )

    hvgs_file = sync_dir / "hvgs.txt"
    if sync_with_other and hvgs_file.exists():
        with hvgs_file.open() as f:
            hvgs = list(map(str.strip, f.readlines()))
            adata_tissue.var["highly_variable"] = False

            # Depends on the df index the user chose
            if adata_tissue.uns["config"]["var_names"] == "gene_ids":
                adata_tissue.var.loc[hvgs, "highly_variable"] = True
            else:
                adata_tissue.var["highly_variable"] = adata_tissue.var["gene_ids"].isin(
                    hvgs
                )
    else:
        # we use n = 2000 since the vignette uses that number
        # We can use flavor as: seurat, seurat_v3 (uses counts), cell_ranger
        # seurat and cell_ranger use log_counts
        sc.pp.highly_variable_genes(adata_tissue, flavor="seurat", n_top_genes=2000)
        if sync_with_other:
            hvgs = adata_tissue.var[adata_tissue.var["highly_variable"]].index.tolist()
            txt_checkpoint(hvgs, filename=sync_dir / "hvgs.txt", include_counter=False)
        # TODO: Add checkpoint?

    # Dimension reduction and clustering

    # This is required since the vignette scales prior to PCA (without actually storing the scaled version)

    # sc.pp.scale(adata_tissue, zero_center=True)
    adata_tissue.X = vp.utils.scale(adata_tissue.X, center=True)
    sc.tl.pca(adata_tissue, use_highly_variable=True, n_comps=30, random_state=1337)

    # TODO: add checkpoint comparing X_pca with reduced_dim(sfe, "PCA")
    mtx_checkpoint(adata_tissue.obsm["X_pca"], filename=checkpoint_dir / "pca.mtx")

    sc.pp.neighbors(adata_tissue, n_pcs=3, use_rep="X_pca")
    sc.tl.leiden(
        adata_tissue,
        random_state=29,
        resolution=0.3,  # to get the same amount of clusters
        key_added="cluster",
    )

    cluster_file = sync_dir / "cluster.txt"

    if sync_with_other:
        if cluster_file.exists():
            pass
        else:
            txt_checkpoint(
                adata_tissue.obs["cluster"].items(),
                filename=cluster_file,
                include_counter=False,
            )

    # TODO: Do we need to compute the umap?
    # sc.tl.umap(adata_tissue)


def download_data(root_dir: Path, data_dir: Path, force: bool = False) -> Path:
    # Download the gene count matrix
    data_dir.mkdir(parents=True, exist_ok=True)

    tar_path_ob = root_dir / "visium_ob.tar.gz"
    url_reads = "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_raw_feature_bc_matrix.tar.gz"
    if not tar_path_ob.exists() or force:
        res = requests.get(url_reads)
        with tar_path_ob.open("wb") as f:
            f.write(res.content)

    # Download the spatial information
    tar_path_sp = root_dir / "visium_ob_spatial.tar.gz"
    url_spatial = "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz"
    if not tar_path_sp.exists() or force:
        res = requests.get(url_spatial)
        with tar_path_sp.open("wb") as f:
            f.write(res.content)

    bc_filenames = ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
    bc_files = [data_dir / "raw_feature_bc_matrix" / fn for fn in bc_filenames]

    spatial_filenames = (
        "aligned_fiducials.jpg",
        "detected_tissue_image.jpg",
        "scalefactors_json.json",
        "spatial_enrichment.csv",
        "tissue_hires_image.png",
        "tissue_lowres_image.png",
        "tissue_positions.csv",
    )

    spatial_files = [data_dir / "spatial" / fn for fn in spatial_filenames]

    if not all(map(Path.exists, bc_files)) or force:
        os.system(f"tar -xvf {tar_path_ob} -C {data_dir}")

    if not all(map(Path.exists, spatial_files)) or force:
        os.system(f"tar -xvf {tar_path_sp} -C {data_dir}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(__name__)
    parser.add_argument("--checkpoint", type=int, choices=(0, 1), default=1)
    parser.add_argument("--datadir", type=str, default="data")
    args = parser.parse_args()

    cwd = Path(__file__).parent
    checkpoint_dir = cwd / f"checkpoints_{args.checkpoint}"
    data_dir = cwd / args.datadir
    print(cwd)
    print(checkpoint_dir)

    main(cwd, checkpoint_dir, data_dir)
