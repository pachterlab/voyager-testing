#!/bin/bash
DIR="."

if [ -n $@ ]; then
	DIR="$1"
fi

mkdir -p ${DIR}/outs
OUTS_DIR="${DIR}/outs"

# Download raw feature matrix
if [ ! -d "${OUTS_DIR}/raw_feature_bc_matrix" ]; then
	COUNTS_URL="https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_raw_feature_bc_matrix.tar.gz"
	curl -s -L "${COUNTS_URL}" | tar xvf - -C "${OUTS_DIR}"
fi

# Download spatial data
if [ ! -d "${OUTS_DIR}/spatial" ]; then
	SPATIAL_URL="https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz"
	curl -s -L "${SPATIAL_URL}" | tar xvf - -C "${OUTS_DIR}"
fi
