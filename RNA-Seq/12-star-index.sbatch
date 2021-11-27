#!/bin/sh -e

##########################################################################
#   https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html
##########################################################################

#SBATCH --cpus-per-task=6
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem=24g
#SBATCH --output=Logs/12-star-index/slurm-%A.out
#SBATCH --error=Logs/12-star-index/slurm-%A.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}
: ${SLURM_CPUS_PER_TASK:=6}

# Document software versions used for publication
uname -a > Logs/12-star-index/os-version-$SLURM_JOB_ID.txt 2>&1
STAR --version > Logs/12-star-index/star-version-$SLURM_JOB_ID.txt 2>&1
samtools --version > Logs/12-star-index/samtools-version-$SLURM_JOB_ID.txt 2>&1
hostname

genome=$(Reference/genome-filename.sh)
gtf=$(Reference/gtf-filename.sh)
printf "Using reference $genome, $gtf...\n"

# Overhang should be read length - 1
STAR --runThreadN $SLURM_CPUS_PER_TASK \
    --runMode genomeGenerate \
    --genomeDir Data/12-star-index \
    --genomeFastaFiles Data/03-reference/$genome \
    --sjdbGTFfile Data/03-reference/$gtf
    --sjdbOverhang 97

ls Data/12-star-index
