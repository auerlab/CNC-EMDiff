#!/bin/sh -e

#
#   Dependencies:
#       Requires SAM FastQC results.  Run after *-qc-sam.sbatch.

if which srun > /dev/null; then
    srun=srun
else
    srun=''
fi

# multiqc: LC_ALL and LANG must be set to a UTF-8 character set
# in your environment in order for the click module to function.
export LC_ALL=en_US.UTF-8

cd Data/11-multiqc-sam
rm -rf *
$srun multiqc --version > ../../Logs/11-multiqc-sam/multiqc-version.txt 2>&1
$srun multiqc ../10-qc-sam 2>&1 | tee ../../Logs/11-multiqc-sam/multiqc.out
