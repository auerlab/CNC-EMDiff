#!/bin/sh -e

if which srun > /dev/null; then
    srun=srun
else
    srun=''
fi

# multiqc: LC_ALL and LANG must be set to a UTF-8 character set
# in your environment in order for the click module to function.
export LC_ALL=en_US.UTF-8

$srun multiqc --version > Logs/11-multiqc-sam/multiqc-version.txt 2>&1

dir=Results/11-multiqc-sam
mkdir -p $dir/Raw $dir/Rmdup

cd $dir/Raw
$srun multiqc ../../10-qc-sam/Raw

cd ../Rmdup
$srun multiqc ../../10-qc-sam/Rmdup
