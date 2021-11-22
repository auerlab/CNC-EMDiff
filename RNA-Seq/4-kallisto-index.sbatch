#!/bin/sh -ex

#SBATCH --mem=4g
#SBATCH --output=Logs/4-kallisto-index/slurm-%A.out
#SBATCH --error=Logs/4-kallisto-index/slurm-%A.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}

# Document software versions used for publication
uname -a > Data/4-kallisto-index/os-version-$SLURM_JOB_ID.txt 2>&1
kallisto version > Data/4-kallisto-index/kallisto-version-$SLURM_JOB_ID.txt 2>&1
samtools --version > Data/4-kallisto-index/samtools-version-$SLURM_JOB_ID.txt 2>&1

reference=$(Reference/reference-filename.sh)
printf "Using reference $reference...\n"

# Needed for kallisto --genomebam
samtools faidx Data/3-reference/$reference

kallisto index --index=Data/4-kallisto-index/all-but-xy.index \
    Data/3-reference/$reference
