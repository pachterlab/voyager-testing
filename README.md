# voyager-testing


This repo tests the implementations of voyager and voyagerpy. To setup testing for a vignette-notebook pair you must take the following steps:

1. For R vignettes, extract the non-plotting code into a script for the corresponding pair
2. Do the same thing for the python notebook, it must reside in the same directory as the R scripts
3. GitHub Actions should take care of the rest (not for the time being though)

The proposed directory structure is as follows:

```
voyager-testing
└── vignettes
    ├── visium_10x
    │   ├── checkpoints_0  # checkpoints for the first run script (e.g. the R script)
    │   │   ├── chkpt_0.mtx
    │   │   ├── ...
    │   │   └── chkpt_n.mtx
    │   ├── checkpoints_1  # checkpoints for the second run script (e.g. the Python script)
    │   │   ├── chkpt_0.mtx
    │   │   ├── ...
    │   │   └── chkpt_n.mtx
    │   ├── download.sh  # bash script to download the data (and unzip).
    │   ├── outs
    │   │   ├── ...
    │   │   └── ...
    │   ├── run.py  # The python script
    │   ├── run.R   # The R script
    │   ├── sync_0  # "Checkpoints" for diverging data from the first run.
    │   │   ├── clusters_0.txt
    │   │   ├── ...
    │   │   └── hvgs_n.txt
    │   └── sync_1  # "Checkpoints" for diverging data, from the second run.
    │       ├── clusters_0.txt
    │       ├── ...
    │       └── hvgs_n.txt
    ├── nonspatial
    │  ├── checkpoints_0/
    │  ├── checkpoints_1/
    │  ├── download.sh
    │  ├── outs/
    │  ├── run.py
    │  ├── run.R
    │  ├── sync_0/
    │  └── sync_1/
    ├── ...
    └── ...
```

To elaborate on the structure above, each vignette-notebook pair has their own directory directly under `vignettes/`.
The R script and Python scripts are directly under this directory. Along with the scripts, there are two `checkpoint_[01]/` directories containing the checkpoint files from each script,
an `outs/` directory which contains the input data for the scripts, and a pair of `sync_[01]/` directory which contains any necessary files to keep the scripts in sync (e.g. gene_var.txt).

The script run first is responsible for writing files in the `sync_0/` and `checkpoints_0/` directories.
The script run second reads from the `sync_0/` directory (if necessary) and writes checkpoints to `checkpoints_1/` directory. It also writes to `sync_1/` directory so we can compare the differences without asserting equality.

The two scripts must write the same checkpoints with the same name (e.g. `ckpt_0.mtx`, `text_file_1.txt`, `pca_2.mtx`, etc). Special matrices which need special consideration, like PCA, must have a name associated with their type, such as `pca_*.mtx`. Further, when comparing set-like data, the name should contain the substring `"set"`. Note, that the numberings of the filenames come automatically. Adding the checkpoints `"ckpt.mtx"`, `"text_file.txt"`, and `"pca.mtx"` (in this order) will create files with the names `"ckpt_0.mtx"`, `"text_file_1.txt"`, `"pca_2.mtx"`. Note that there is a separate counter for the directories. Thus, if the `"text_file.txt"` were to be added as a sync point the names would be `"ckpt_0.mtx"`, `"pca_1.mtx"` in the checkpoints directory, and `"text_file_0.txt"` in the sync directory.


## Setting up a compatibility test of a new vignette

In order to compare the python an R vignettes, these steps must be taken:

1. Create a unique directory name in `vignettes/`
2. Add the name of the new directory to `manifest.txt`
3. Convert the R vignette to an executable `run.R` script in the new directory.
4. Convert the Python notebook to an executable `run.py` script in the new directory.
5. (Optional) Create an executable `download.sh` bash script in the new directory. The download script will be called with the directory it lives in. This script should be written such that it downloads and decompresses the input data into a directory (e.g. `outs/`) residing in this vignette dir.
6. Run `./make_tests.sh` in the root directory of this project. You can comment out (with `#`)other tests in the manifest file when running locally.

### Converting vignettes to `run.R`

When converting a vignette to an R script, there are a few steps:

1. Copy the runnable code to the empty run.R file
2. Remove all non-essential code, e.g. plotting and printing/displaying, or any computation that ONLY affect plots.
3. To prevent messages from the import statements, wrap the library imports in `suppressMessages()`.
4. Add the following lines (you can tweak them however you want) after the imports: 
```R
root_path <- normalizePath(Sys.getenv("RPATH"))
utils_path <- file.path(root_path, "utils.R")
source(utils_path)
init.dirs()
```
5. Adjust the file path used for reading the input data.
6. (Optional) Remove the code responsible for downloading and preparing the data, if you have provided a `download.sh` script.
7. Add the checkpoints you need to check.

The `utils.R` script provides the following functions and variables:

- `init.dirs()`: reads from the command line and sets the variables listed below.
- `normalize.filename(filename, ...)`: returns the path to the relevant sync/checkpoint directory, with the relevant counter, and increases said counter by one. This is called in each `checkpoint.*` function.
- `checkpoint.csv(filename, dat, ...)`: Write a csv checkpoint/sync. dat must be castable as data.frame.
- `checkpoint.txt(filename, dat, ...)`: Write a txt checkpoint/sync, defaulting to a single column.
- `checkpoint.mtx(filename, dat, ...)`: Write a sparse matrix to an mtx checkpoint/sync.
- `checkpoint.knn(filename, dat, ...)`: Converts a `listw` object to sparse matrix and writes an mtx checkpoint/sync.
- `.checkpoint.counter`, `.sync.counter`: used and updated in `normalize.filename()` for counting the number of files in the checkpoints/sync directories. Should not be touched unless you know what you're doing.
- `.is.primary`: TRUE if this is the primary (e.g. first) run between the two runs. This is set in `init.dirs()`.
- `.root.dir`: the directory under for this test.
- `.checkpoint.dir`: The checkpoint directory
- `.sync.dir`: The sync directory for this run.

To write a `.csv` checkpoint of a dataframe `df`:

```R
checkpoint.csv("filename.csv", df) # sync=FALSE is default
```

To write a `.csv` sync point of the same data frame:
```R
checkpoint.csv("filename.csv", df, sync=TRUE)
```

**NOTE:** You should not insert the counter yourself, this is done automatically when each of the writing functions call `normalize.filename(filename)`.

### Convert an ipynb notebook to `run.py`

**Tip**: For a cleaner workspace you can start with step 2 and convert the notebook as the last step. 

The steps are somewhat similar to the process above:

1. Copy the runnable code from the ipynb. You can also export the notebook using the interface of `jupyter`, or through the command line `jupyter nbconvert notebook.ipynb --to python`.
2. Remove all non-essential code e.g. plotting, printing/displaying.
3. Add the following lines after the imports:
```python
from utils import Checkpoint
checkpoint = Checkpoint()
# Set the relevant directory. This depends on the location of the downloaded data.
outs_dir = pathlib.Path(checkpoint.root_dir / "outs")
# Use this directory as input to the VoyagerPy read* function.
```
4. Adjust the file fath used for reading the input data.
5. (Optional) Remove the code responsible for downloading and preparing the data, if you have provided the `download.sh` script.
6. Add the checkpoints you need to check.

The `utils.py` module provides the `Checkpoint` class. Once initialized, you can add a checkpoint via:

```python
# checkpoint = Checkpoint()
checkpoint.add("filename.csv", df, type="csv")
```

To add a sync point:

```python
loaded_df = checkpoint.sync("filename.csv", df, "csv")
```

**NOTE:** You should only give the name of the file sans numbering (the checkpoint class handles that for you). Numbering doesn't affect it's behaviour though.


