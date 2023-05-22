import numpy as np
import os
import pandas as pd
from collections import defaultdict
from pathlib import Path
from scipy import io, sparse
from typing import Any, Union, Optional, Sequence, List, Literal
import argparse

FileType = Literal["txt", "mtx", "csv"]


class Checkpoint(object):
    def __init__(self, argstr=None) -> None:
        args = self._parse_args(argstr)
        self.primary = bool(args.primary)

        self.root_dir = Path(args.dir).resolve()
        self.sync_dir = self._get_sync_path(self.primary)
        self.checkpoint_dir = self.root_dir / f"checkpoints_{int(not self.primary)}"

        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self.sync_dir.mkdir(parents=True, exist_ok=True)

        self.quiet = args.quiet

        self.counter = defaultdict(int)

    def _get_sync_path(self, primary: bool):
        return self.root_dir / f"sync_{int(not primary)}"

    def _parse_args(self, argstr=None):
        parser = argparse.ArgumentParser()
        parser.add_argument("-p", "--primary", action="store_true")
        parser.add_argument("-d", "--dir", type=str, required=True)
        parser.add_argument("-q", "--quiet", action="store_true")
        args = parser.parse_known_args(argstr)[0]
        return args

    def sync(
        self,
        filename: str,
        data: Any = None,
        file_type: Optional[FileType] = None,
        write_kwargs: Optional[dict] = None,
        **kwargs,
    ):
        filepath = self.normalize_path(filename, sync=True)

        if file_type is None:
            file_type = filepath.suffix.strip(".")  # type: ignore

        if data is None:
            raise ValueError("Must provide data to any sync point.")

        self._write(filepath, data, file_type, **(write_kwargs or {}))

        if not self.primary:
            read_path = self._get_sync_path(not self.primary) / filepath.name

            # load data
            if file_type == "mtx":
                return self.read_mtx(read_path, **kwargs)
            elif file_type == "txt":
                return self.read_txt(read_path, **kwargs)
            elif file_type == "csv":
                return self.read_csv(read_path, **kwargs)
            else:
                raise ValueError(f"Invalid file type '{file_type}'")

        return data

    def normalize_path(self, filename, sync: bool = False):
        if sync:
            file_path = self.sync_dir / filename
            i_file = self.counter[self.sync_dir]
            self.counter[self.sync_dir] += 1
        else:
            file_path = self.checkpoint_dir / filename
            i_file = self.counter[self.checkpoint_dir]
            self.counter[self.checkpoint_dir] += 1

        suffixes = "".join(file_path.suffixes)
        file_path = (
            file_path.with_name(f"{file_path.stem}_{i_file}")
            .with_suffix(suffixes)
            .resolve()
        )
        if not self.quiet:
            print(file_path.relative_to(self.root_dir))

        return file_path

    def _write(self, filepath: Path, data: Any, type: Optional[str] = None, **kwargs):
        if type is None:
            type = filepath.suffix.strip(".")

        filepath = filepath.with_suffix(f".{type}")
        if type == "mtx":
            self.write_mtx(filepath, data, **kwargs)
        elif type == "txt":
            self.write_txt(filepath, data, **kwargs)
        elif type == "csv":
            self.write_csv(filepath, data, **kwargs)
        else:
            raise ValueError(f"Invalid type '{type}'")

    def add(
        self,
        filename: str,
        data: Any,
        type: Optional[FileType] = None,
        sync: bool = False,
        **kwargs,
    ):
        """Write a checkpoint file

        Parameters
        ----------
        filename : str
            The name of the file to write
        data : Any
            The data to write to the file
        type : Optional[str], optional
            The type of file to write, by default None

        Returns
        -------
        None
        """
        if sync:
            return self.sync(filename, data, type, **kwargs)

        filepath = self.normalize_path(filename)

        self._write(filepath, data, type, **kwargs)

    def read_mtx(self, filename: Path, **kwargs):
        return io.mmread(filename)

    def read_txt(self, filename: Path, **kwargs):
        with filename.open("rt") as f:
            return [line.strip() for line in f]

    def read_csv(self, filename: Path, **kwargs) -> pd.DataFrame:
        return pd.read_csv(filename, **kwargs)

    def write_csv(self, filename: Path, df, **kwargs):
        if isinstance(df, pd.DataFrame):
            df.to_csv(filename, **kwargs)
        elif isinstance(df, np.ndarray):
            pd.DataFrame(df).to_csv(filename, **kwargs)
        else:
            raise ValueError(f"Invalid type '{type(df)}' for df")

    def write_mtx(self, filename: Path, df, **kwargs):
        io.mmwrite(filename, df)

    def write_txt(
        self,
        filename: Path,
        items: List[Union[str, Sequence[str]]],
        sep: str = ",",
        **kwargs,
    ):
        with filename.open("wt") as f:
            rows = items if isinstance(items[0], str) else map(sep.join, items)
            f.writelines([f"{row}\n" for row in rows])
