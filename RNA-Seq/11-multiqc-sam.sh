#!/bin/sh -e

if which srun; then
    srun=srun
else
    srun=''
fi

$srun multiqc --version > Data/10-qc-sam/multiqc-version.txt 2>&1

cd Data/10-qc-sam && $srun multiqc .
