#!/bin/sh -e

#SBATCH --mem=1g
#SBATCH --output=Logs/08-hisat2-index/slurm-%A.out
#SBATCH --error=Logs/08-hisat2-index/slurm-%A.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}

# Document software versions used for publication
uname -a > Logs/08-hisat2-index/os-version-$SLURM_JOB_ID.txt 2>&1
hisat2 --version > Logs/08-hisat2-index/hisat2-version-$SLURM_JOB_ID.txt 2>&1
samtools --version > Logs/08-hisat2-index/samtools-version-$SLURM_JOB_ID.txt 2>&1
hostname

# Run hisat2-build on a copy in 08-hisat2-index so it will put the .ht2
# files there
genome=$(Reference/genome-filename.sh)
ln -f Data/03-reference/$genome Data/08-hisat2-index
genome=Data/08-hisat2-index/$genome
printf "Using reference $genome...\n"

if [ ! -e $genome.8.ht2 ]; then
    printf "Building $genome.*.ht2...\n"
    hisat2-build $genome $genome
fi
if [ ! -e $genome.fai ]; then
    printf "Building $genome.fai...\n"
    samtools faidx $genome
fi
ls Data/08-hisat2-index
