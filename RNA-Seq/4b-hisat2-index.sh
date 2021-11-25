#!/bin/sh -e

#SBATCH --mem=1g
#SBATCH --output=Logs/4b-hisat2-index/slurm-%A.out
#SBATCH --error=Logs/4b-hisat2-index/slurm-%A.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}

# Document software versions used for publication
uname -a > Logs/4b-hisat2-index/os-version-$SLURM_JOB_ID.txt 2>&1
hisat2 --version > Logs/4b-hisat2-index/hisat2-version-$SLURM_JOB_ID.txt 2>&1
samtools --version > Logs/4b-hisat2-index/samtools-version-$SLURM_JOB_ID.txt 2>&1

reference=Data/3-reference/$(Reference/reference-filename.sh)
cp $reference Data/4b-hisat2-index
reference=Data/4b-hisat2-index/$(basename $reference)
printf "Using reference $reference...\n"

if [ ! -e $reference.8.ht2 ]; then
    printf "Building $reference.*.ht2...\n"
    hisat2-build $reference $reference
fi
if [ ! -e $reference.fai ]; then
    printf "Building $reference.fai...\n"
    samtools faidx $reference
fi
ls Data/4b-hisat2-index
