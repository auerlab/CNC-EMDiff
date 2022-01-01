#!/bin/sh -e

if which srun; then
    srun=srun
else
    srun=''
fi

$srun multiqc --version > Data/02-qc-raw/multiqc-version.txt 2>&1

cd Data/02-qc-raw && $srun multiqc .
