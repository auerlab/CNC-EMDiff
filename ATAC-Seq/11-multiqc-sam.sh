#!/bin/sh -e

if which srun; then
    srun=srun
else
    srun=''
fi

$srun multiqc --version > 5-qc-sam/multiqc-version.txt 2>&1

(cd 5-qc-sam/Raw && $srun multiqc .)
(cd 5-qc-sam/Trimmed && $srun multiqc .)
