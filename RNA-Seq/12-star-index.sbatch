#!/bin/sh -e

##########################################################################
#   https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html
##########################################################################

#SBATCH --mem=21g
#SBATCH --ntasks=6
#SBATCH --output=Logs/12-star-index/slurm-%A.out
#SBATCH --error=Logs/12-star-index/slurm-%A.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}
: ${SLURM_NTASKS:=6}

# Document software versions used for publication
uname -a > Logs/12-star-index/os-version-$SLURM_JOB_ID.txt 2>&1
STAR --version > Logs/12-star-index/star-version-$SLURM_JOB_ID.txt 2>&1
samtools --version > Logs/12-star-index/samtools-version-$SLURM_JOB_ID.txt 2>&1

genome=$(Reference/genome-filename.sh)
gtf=$(Reference/gtf-filename.sh)
printf "Using reference $genome, $gtf...\n"

# Overhang should be read length - 1
STAR --runThreadN $SLURM_NTASKS \
    --runMode genomeGenerate \
    --genomeDir Data/12-star-index \
    --genomeFastaFiles Data/03-reference/$genome \
    --sjdbGTFfile Data/03-reference/$gtf
    --sjdbOverhang 97

ls Data/12-star-index
