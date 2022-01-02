#!/bin/sh -e

if which srun > /dev/null; then
    srun=srun
else
    srun=''
fi

# multiqc: LC_ALL and LANG must be set to a UTF-8 character set
# in your environment in order for the click module to function.
export LC_ALL=en_US.UTF-8

cd Data/06-multiqc-trimmed
rm -rf *
$srun multiqc --version > multiqc-version.txt 2>&1
$srun multiqc ../05-qc-trimmed
