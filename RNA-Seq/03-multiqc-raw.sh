#!/bin/sh -e

#
#   Dependencies:
#       Requires raw FastQC results.  Run after *-qc-raw.sbatch.

if which srun > /dev/null; then
    srun=srun
else
    srun=''
fi

# multiqc: LC_ALL and LANG must be set to a UTF-8 character set
# in your environment in order for the click module to function.
export LC_ALL=en_US.UTF-8

cd Results/04-multiqc-raw
rm -rf *
$srun multiqc --version > ../../Logs/04-multiqc-raw/multiqc-version.txt 2>&1
$srun multiqc ../02-qc-raw 2>&1 | tee ../../Logs/04-multiqc-raw/multiqc.out
