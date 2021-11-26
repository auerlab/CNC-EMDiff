#!/bin/sh -e

#SBATCH --mem=4g
#SBATCH --output=Logs/04-kallisto-index/slurm-%A.out
#SBATCH --error=Logs/04-kallisto-index/slurm-%A.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}

# Document software versions used for publication
uname -a > Logs/04-kallisto-index/os-version-$SLURM_JOB_ID.txt 2>&1
kallisto version > Logs/04-kallisto-index/kallisto-version-$SLURM_JOB_ID.txt 2>&1
samtools --version > Logs/04-kallisto-index/samtools-version-$SLURM_JOB_ID.txt 2>&1

transcriptome=Data/03-reference/$(Reference/transcriptome-filename.sh)
printf "Using reference $transcriptome...\n"

# Needed for kallisto --genomebam
if [ ! -e $transcriptome.fai ]; then
    printf "Building $transcriptome...\n"
    samtools faidx $transcriptome
fi

printf "Building kallisto index...\n"
set -x
kallisto index --index=Data/04-kallisto-index/all-but-xy.index $transcriptome
