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
    │   ├── data
    │   │   ├── ...
    │   │   └── ...
    │   ├── notebook.py  # The python script
    │   ├── vignette.R  # The R script
    │   └── sync_data  # "Checkpoints" for diverging data, if the second run scripts needs to read these
    │       ├── clusters.txt
    │       └── hvgs.txt
    ├── visium_spatial
    │  ├── ...
    │  └── ...
    ├── ...
    └── ...
```

To elaborate on the structure above, each vignette-notebook pair has their own directory directly under `vignettes`.
The R script and Python scripts are directly under this directory. Along with the scripts, there are two `checkpoint` directories containing the checkpoint files from each script,
a `data` directory which contains the input data for the scripts, and an optional `sync_data` directory which contains any necessary files to keep the scripts in sync (e.g. hvgs.txt).

The script run first is responsible for writing files in the `sync_data` directory, as well as the `checkpoints_0` directory.
The script run second reads from the `sync_data` directory (if necessary) and writes checkpoints to `checkpoints_1` directory.

The two scripts must write the same checkpoints with the same name (or at least the same numbering e.g. `ckpt_0.mtx`, `text_file_1.txt`, `pca_2.mtx`, etc). Special matrices which need special consideration, like PCA, must have a name associated with their type, such as `pca_*.mtx`.