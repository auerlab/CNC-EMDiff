#!/bin/sh -e

#
#   Dependencies:
#       Requires hisat2 alignments.  Run after *-hisat2-align.sbatch.

# Computing abundances takes a little while, so we do this in parallel
# One core for fasda abundance, one for samtools
#SBATCH --array=1-18
#SBATCH --cpus-per-task=2
#SBATCH --output=Logs/19-fasda-abundance-hisat2/slurm-%A_%a.out
#SBATCH --error=Logs/19-fasda-abundance-hisat2/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
fasda --version
samtools --version
pwd

hisat2_dir=Results/17-hisat2-align
fasda_dir=Results/19-fasda-abundance-hisat2

# Need to fetch gff3 for computing abundances
gff_filename=$(Reference/gff-filename.sh)
reference_dir=Results/07-reference

##########################################################################
#   Compute abundances
##########################################################################

file=`ls $hisat2_dir/*-sample$SLURM_ARRAY_TASK_ID-*.bam`

fasda abundance \
    --output-dir $fasda_dir \
    $reference_dir/$gff_filename \
    $file

ab=${file%.bam}-abundance.tsv
head $ab
