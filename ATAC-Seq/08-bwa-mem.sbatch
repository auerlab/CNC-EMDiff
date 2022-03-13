#!/bin/sh -e

##########################################################################
#   Script description:
#       Align reads to reference genome
#
#   Usage:
#       SLURM cluster:
#           sbatch 08-bwa-mem.sbatch
#
#   History:
#       Based on the work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#   Date        Name        Modification
#   2020-02-20  Jason Bacon Begin
##########################################################################

##########################################################################
#   Main
##########################################################################

# SLURM parameters
# Not yet sure how many threads are useful.  Maybe use --exclusive and all
# available cores? Smaller thread counts are usually more efficient as
# they will require less communication overhead.  If cores are limited,
# this will get the job done faster.
#SBATCH --array=1-18 --cpus-per-task=4
# Based on 2 threads on Albacore cluster
# 12 threads on Peregrine used close to 6g
#SBATCH --mem=6g
#SBATCH --output=Logs/08-bwa-mem/bwa-mem-%A_%a.out
#SBATCH --error=Logs/08-bwa-mem/bwa-mem-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_CPUS_PER_TASK:=4}

# Document software versions used for publication
uname -a
bwa || true
pwd

ref_file=$(../RNA-Seq/Reference/genome-filename.sh)

# One iteration if running under SLURM, all iterations otherwise
cd Data/08-bwa-mem
printf "Running task $SLURM_ARRAY_TASK_ID...\n"
infile1=$(ls ../04-trim/*-sample$SLURM_ARRAY_TASK_ID-*-R1.fastq.gz)
infile2=$(ls ../04-trim/*-sample$SLURM_ARRAY_TASK_ID-*-R2.fastq.gz)
outfile=$(basename ${infile1%-R1.fastq.gz}).sam

set -x
bwa mem -M -t $SLURM_CPUS_PER_TASK \
    ../07-bwa-index/$ref_file $infile1 $infile2 > $outfile
