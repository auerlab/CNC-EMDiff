#!/bin/sh -e

##########################################################################
#   Script description:
#       QC BWA alignments
#
#   Usage:
#       SLURM cluster:
#           sbatch 10-qc.sbatch
#
#   History:
#       Based on the work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#   Date        Name        Modification
#   2020-02-25  Jason Bacon Begin
##########################################################################

##########################################################################
#   Main
##########################################################################

#SBATCH --array=1-18
# FIXME: Check actual memory use
#SBATCH --mem=1g
#SBATCH --output=Logs/10-qc-sam/slurm-%A_%a.out
#SBATCH --error=Logs/10-qc-sam/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
fastqc --version
pwd

qc_dir=Data/10-qc-sam
mkdir -p $qc_dir/Raw
infile=Data/08-bwa-mem/*-sample${SLURM_ARRAY_TASK_ID}-*.sam
printf "Processing $infile...\n"
fastqc $infile -o $qc_dir/Raw

mkdir -p $qc_dir/Rmdup
infile=Data/09-remove-duplicates/*-sample${SLURM_ARRAY_TASK_ID}-*-mapq1.bam
printf "Processing $infile...\n"
fastqc $infile -o $qc_dir/Rmdup
