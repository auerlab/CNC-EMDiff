#!/bin/sh -e

##########################################################################
#   Description:
#       Run STAR aligner on each RNA sample.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#       
#   History:
#   Date        Name        Modification
#   2021-11-26  Jason Bacon Begin
##########################################################################

##########################################################################
#   https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html
##########################################################################

#SBATCH --array=1-18
# 1 for STAR, 1 for each xzcat
#SBATCH --cpus-per-task=3
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem=24000
#SBATCH --output=Logs/13-star-align/slurm-%A_%a.out
#SBATCH --error=Logs/13-star-align/slurm-%A_%a.err

##############################################################################
# Align with star, which can handle splice junctions in RNA reads

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}
: ${SLURM_CPUS_PER_TASK:=4}

# Document software versions used for publication
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
    uname -a > Logs/13-star-align/os-version-$SLURM_JOB_ID.txt 2>&1
    STAR --version > Logs/13-star-align/star-version-$SLURM_JOB_ID.txt 2>&1
fi
hostname

build=$(../Common/genome-build.sh)
release=$(../Common/genome-release.sh)
genome=$(Reference/genome-filename.sh)

xz1=$(echo Data/01-trim/*-sample${SLURM_ARRAY_TASK_ID}-*-R1.fastq.xz)
xz2=$(echo Data/01-trim/*-sample${SLURM_ARRAY_TASK_ID}-*-R2.fastq.xz)

xzb=$(basename $xz1)
filename_prefix=${xzb%-R*}-
bam=Data/13-star-align/${filename_prefix}Aligned.sortedByCoord.out.bam

set -x
STAR --genomeDir Data/12-star-index \
    --runThreadN 1 \
    --readFilesIn $xz1 $xz2 \
    --readFilesCommand xzcat \
    --outFileNamePrefix Data/13-star-align/$filename_prefix \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard

samtools index $bam
