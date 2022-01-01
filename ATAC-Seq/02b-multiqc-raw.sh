#!/bin/sh -e

if which srun; then
    srun=srun
else
    srun=''
fi

$srun multiqc --version > 2-qc/multiqc-version.txt 2>&1

(cd 2-qc/Raw && $srun multiqc .)
(cd 2-qc/Trimmed && $srun multiqc .)
