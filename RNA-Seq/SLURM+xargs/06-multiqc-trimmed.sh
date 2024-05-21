#!/bin/sh -e

#
#   Dependencies:
#       Requires trimmed FastQC results.  Run after *-qc-trimmed.sbatch.

if which srun > /dev/null; then
    srun=srun
else
    srun=''
fi

# multiqc: LC_ALL and LANG must be set to a UTF-8 character set
# in your environment in order for the click module to function.
export LC_ALL=en_US.UTF-8

cd Results/06-multiqc-trimmed
rm -rf *
$srun multiqc --version > ../../Logs/06-multiqc-trimmed/multiqc-version.txt 2>&1
$srun multiqc ../05-qc-trimmed 2>&1 | tee ../../Logs/06-multiqc-trimmed/multiqc.out
