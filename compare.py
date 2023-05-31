"""Compare sync points and checkpoints between implementations of the same vignette/notebook.
"""

import argparse
import filecmp
import pandas as pd

import numpy as np
from scipy import io, sparse as sp
from pathlib import Path
from typing import Tuple, Callable


def compare_txt(file_0: Path, file_1: Path, verbose: bool = False, **kwargs) -> bool:
    """Compare two text files line by line.

    Args:
        file_0: Path to first file.
        file_1: Path to second file.
        verbose: Print differences if True.

    Returns:
        True if files are identical, False otherwise.
    """
    with file_0.open("r") as f0, file_1.open("r") as f1:
        lines_0 = list(f0)
        lines_1 = list(f1)
        if len(lines_0) != len(lines_1):
            return False
        for line0, line1 in zip(lines_0, lines_1):
            if line0 != line1:
                return False
    return True


# function to compare the set differences of two files
def compare_set(file_0: Path, file_1: Path, verbose: bool = False, **kwargs) -> bool:
    """"""
    with file_0.open("r") as f0, file_1.open("r") as f1:
        lines_0 = set(f0)
        lines_1 = set(f1)
        if lines_0 != lines_1 and verbose:
            n_mismatch = len(lines_0 - lines_1)
            max_lines = max(verbose * 10, n_mismatch)
            print(f"{file_0.parent.name} - {file_1.parent.name}:")
            print("", *sorted(lines_0 - lines_1)[:max_lines], sep="    ")

            if max_lines < n_mismatch:
                print(f"    {n_mismatch - max_lines} more lines")
            print(f"    (total: {n_mismatch})")

            n_mismatch = len(lines_1 - lines_0)
            max_lines = max(verbose * 10, n_mismatch)
            print(f"{file_1.parent.name} - {file_0.parent.name}:")
            print("", *sorted(lines_1 - lines_0)[:max_lines], sep="    ")
            if max_lines < n_mismatch:
                print(f"    {n_mismatch - max_lines} more lines")
            print(f"(total: {n_mismatch})")

            return False
    return True


def compare_mtx(
    file_0: Path, file_1: Path, eps: float = 1e-6, verbose: bool = False, **kwargs
) -> bool:
    mtx_0 = io.mmread(file_0)
    mtx_1 = io.mmread(file_1)

    if sp.issparse(mtx_0) != sp.issparse(mtx_1):
        return False
    if not (mtx_0.shape == mtx_1.shape or mtx_0.T.shape == mtx_1.shape):
        return False

    if mtx_0.shape == mtx_1.T.shape:
        mtx_1 = mtx_1.T
    mtx_0 = mtx_0.A
    mtx_1 = mtx_1.A

    res = np.allclose(mtx_0, mtx_1, atol=eps)

    if verbose:
        print("max_diff: ", np.abs(mtx_0 - mtx_1).max())
    return res


def compare_pca(
    file_0: Path,
    file_1: Path,
    pca_eps: float = 1e-6,
    verbose: bool = False,
    test_ortho: bool = False,
    **kwargs,
) -> bool:
    m1, m2 = io.mmread(file_0), io.mmread(file_1)
    if sp.issparse(m1):
        m1 = m1.A
    if sp.issparse(m2):
        m2 = m2.A

    assert m1.shape == m2.shape
    assert m1.dtype == m2.dtype

    success = True

    max_diff_parallel = 0
    max_diff_orthogonal = 0

    self_max_parallel = 0
    self_max_orthogonal = 0
    max_diff_len = 0.0

    for i in range(m1.shape[1]):
        n1 = np.linalg.norm(m1[:, i])
        n2 = np.linalg.norm(m2[:, i])

        d = np.dot(m1[:, i], m2[:, i]) / (n1 * n2)
        self_d = np.dot(m1[:, i], m1[:, i]) / (n1 * n1)
        max_diff_len = max(max_diff_len, np.abs(n1 - n2))

        if not (1 - pca_eps <= abs(d) <= 1 + pca_eps):
            diff = abs(1 - abs(d))
            if verbose > 1:
                print(diff)
            success = False

        max_diff_parallel = max(max_diff_parallel, abs(1 - abs(d)))
        self_max_parallel = max(self_max_parallel, abs(1 - abs(self_d)))

        for j in range(i):
            n2 = np.linalg.norm(m2[:, j])
            d = np.dot(m1[:, i], m2[:, j]) / (n1 * n2)
            self_d = np.dot(m1[:, i], m1[:, j]) / (n1 * n1)
            if abs(d) > pca_eps and test_ortho:
                success = False

            max_diff_orthogonal = max(max_diff_orthogonal, abs(d))
            self_max_orthogonal = max(self_max_orthogonal, abs(self_d))

    if verbose:
        print()
        ortho_str = "" if test_ortho else ", not tested"
        if verbose > 2:
            print(f"\tmax PCA diff (self parallel): {self_max_parallel}")
            print(f"\n\tmax PCA diff (self orthogonal): {self_max_orthogonal}")

        print(f"\tmax PCA diff (parallel): {max_diff_parallel}")
        print(f"\tmax PCA diff (orthogonal{ortho_str}): {max_diff_orthogonal}")
        print(f"\tmax PCA diff (length): {max_diff_len}")

    return success


def compare_csv(
    file_0: Path, file_1: Path, eps: float = 1e-6, verbose: bool = False, **kwargs
) -> bool:
    df_0 = pd.read_csv(file_0, index_col=0).sort_index()
    df_1 = pd.read_csv(file_1, index_col=0).sort_index()

    if df_0.shape != df_1.shape:
        if verbose:
            print("    shape mismatch!")
        return False

    if df_0.columns.tolist() != df_1.columns.tolist():
        if verbose:
            print("    column name mismatch!")
        return False
    if df_0.index.tolist() != df_1.index.tolist():
        if verbose:
            print("    index name mismatch!")
        return False

    all_close = np.allclose(df_0.values, df_1.values, equal_nan=True, atol=eps)

    if not all_close and verbose:
        comp_df = pd.DataFrame(
            index=df_0.columns, columns=["close", "max_eps", "mean_eps"]
        )
        for col in df_0.columns:
            col_0 = df_0[col]
            col_1 = df_1[col]
            cols_close = np.allclose(
                col_0.values, col_1.values, equal_nan=True, atol=eps
            )
            max_eps = np.nanmax(np.abs(col_0.values - col_1.values))
            mean_eps = np.nanmean(np.abs(col_0.values - col_1.values))
            comp_df.loc[col] = [cols_close, max_eps, mean_eps]

        print(comp_df)
    return all_close


def get_common_files(dir_0, dir_1, verbose: int = 0):
    res = filecmp.dircmp(dir_0, dir_1)
    if res.left_only and verbose > 0:
        print(f"Files not found in {dir_1}:", *res.left_only, sep="\n\t")
    if res.right_only and verbose > 0:
        print(f"Files not found in {dir_0}:", *res.right_only, sep="\n\t")

    def comp(fn):
        stem = fn.split(".")[0]
        return int(stem.rsplit("_", 1)[1])

    return sorted(res.common_files, key=comp)


def get_compare_func(filename: str) -> Callable[[Path, Path, bool], bool]:
    if filename.endswith(".csv"):
        return compare_csv
    elif "pca" in filename:
        return compare_pca
    elif filename.endswith(".mtx"):
        return compare_mtx
    elif filename.endswith(".txt"):
        if "_set" in filename:
            return compare_set
        return compare_txt
    else:
        raise ValueError(f"unknown file type: {filename}")


def test_files(dir_0: Path, dir_1: Path, **kwargs) -> Tuple[bool, pd.DataFrame]:
    verbose = kwargs.get("verbose", 1)

    ret = True
    common_files = get_common_files(dir_0, dir_1, verbose=verbose)
    comp_df = pd.DataFrame(index=common_files)
    for filename in common_files:
        file_0 = dir_0 / filename
        file_1 = dir_1 / filename

        if verbose:
            print("Comparing", filename, end=":\n")
        compare = get_compare_func(filename)
        res = compare(file_0, file_1, **kwargs)  # type: ignore
        ret = ret and res
        if verbose:
            print("Equal:", res)
        comp_df.loc[filename, "equal"] = res

    return ret, comp_df


def parse_args():
    parser = argparse.ArgumentParser(
        "compare",
        description=__doc__,
    )
    parser.add_argument(
        "dir",
        type=Path,
        help="Directory where the syncdirs and checkpoints are stored.",
    )
    parser.add_argument(
        "--eps", type=float, default=1e-6, help="Epsilon for comparing matrices."
    )
    parser.add_argument(
        "--pca-eps",
        type=float,
        default=1e-6,
        help="Epsilon for comparing PCA embeddings. Default 1e-6.",
    )
    parser.add_argument(
        "--sync",
        type=int,
        default=1,
        help="Whether to compare the syncdirs (1) or checkpoints (0). Default 1.",
    )

    parser.add_argument(
        "-v", "--verbose", action="count", default=0, help="Verbose output."
    )

    return parser.parse_args()


def main(args):
    sync_dirs = tuple(args.dir / f"sync_{i}" for i in range(2))

    checkpoint_dirs = tuple(args.dir / f"checkpoints_{i}" for i in range(2))
    kwargs = vars(args)
    kwargs.pop("dir")

    sync_df = None
    comp_df = None
    if all(d.exists() for d in sync_dirs) and args.sync:
        res, sync_df = test_files(*sync_dirs, **kwargs)
        print("-----------------------------")
        print("Syncdirs are", "equal" if res else "not equal")
        print()

    if all(d.exists() for d in checkpoint_dirs):
        res, comp_df = test_files(*checkpoint_dirs, **kwargs)
        print("-----------------------------")
        print("Checkpoints are", "equal" if res else "not equal")

        print()

    print("Summary:")
    if sync_df is not None:
        print("Sync")
        print(sync_df)
    if comp_df is not None:
        print("Checkpoints")
        print(comp_df)


if __name__ == "__main__":
    args = parse_args()
    main(args)
