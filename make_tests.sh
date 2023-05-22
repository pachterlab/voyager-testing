#!/bin/bash

set -e
CWD="$(dirname $0)"

export PYTHONPATH="$(dirname $0):$PYTHONPATH"
export RPATH="$(dirname $0)"

PROGRAM_1="Rscript"
PROGRAM_2="python"

FILE_1="run.R"
FILE_2="run.py"

grep -v '^#' ${CWD}/vignettes/manifest.txt | while read -r vignette; do
	DIR="${CWD}/vignettes/${vignette}"
	if [ ! -d "$DIR" ]; then
		echo "Directory ${DIR} not found... skipping"
	else
		# Download data
		if [ -f "$DIR/download.sh" ]; then
			bash $DIR/download.sh $DIR
		fi

		$PROGRAM_1 "${DIR}/${FILE_1}" --primary --dir="${DIR}"
		$PROGRAM_2 $DIR/${FILE_2} --dir="${DIR}"

		python -m compare "${DIR}" --sync=1 -v
	fi

done
