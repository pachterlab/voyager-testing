#!/bin/bash
DIR="."

if [ -n $@ ]; then
	DIR="$1"
fi

OUTS_DIR="${DIR}/outs"

# Download filtered feature matrix
COUNTS_URL="https://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3_nextgem/5k_pbmc_v3_nextgem_filtered_feature_bc_matrix.tar.gz"
mkdir -p "${DIR}/outs"

if [ ! -d "${OUTS_DIR}/filtered_feature_bc_matrix" ]; then
	curl -s -L "${COUNTS_URL}" | tar xvf - -C "${OUTS_DIR}"
fi
