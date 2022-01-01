#!/bin/sh -e

if which srun; then
    srun=srun
else
    srun=''
fi

$srun multiqc --version > Data/04-qc-trimmed/multiqc-version.txt 2>&1

cd Data/04-qc-trimmed && $srun multiqc .
