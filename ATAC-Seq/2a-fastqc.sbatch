#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on raw and trimmed data for comparison
#
#   Usage:
#       SLURM cluster:
#           sbatch 2-qc.sbatch
#       No cluster:
#           ../slurm-sim 2-qc.sbatch |& tee 2.log
#
#   History:
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#   Date        Name        Modification
#   2020-02-05  Jason Bacon Begin
##########################################################################

##########################################################################
#   Main
##########################################################################

# SLURM parameters
# --cpus-per-task is not ideal.  It requires all cores on the same node,
# where we really need 2 pairs.
#SBATCH --array=1-18 --cpus-per-task=4
#SBATCH --mem=1g
#SBATCH --output=2-qc/slurm-%A_%a.out
#SBATCH --error=2-qc/slurm-%A_%a.err

# SLURM will set SLURM_ARRAY_TASK_ID to one of the values in --array above.
# If not running on a SLURM cluster, loop through all the values by spoofing
# SLURM_ARRAY_TASK_ID to the entire list.  Make sure these variables agree
# with the #SBATCH directives above and cutadapt command flags below.
# To run specific samples, set these variables before running the script,
# e.g.
#
#   env SLURM_ARRAY_TASK_ID='1 5' ./2-qc.sbatch

# Default if not set by SLURM or user
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_CPUS_PER_TASK:=2}

# Document software versions used for publication
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
    uname -a > 2-qc/os-version.txt 2>&1
    fastqc --version > 2-qc/fastqc-version.txt 2>&1
fi

mkdir -p 2-qc/Raw 2-qc/Trimmed

# One iteration if running under SLURM, all iterations otherwise
printf "Running task $SLURM_ARRAY_TASK_ID...\n"

# Raw files from sequencing center
raw1=$(ls ../Raw/190822_AHFN3KDRXX-Lane1/*_S${SLURM_ARRAY_TASK_ID}_*R1*.fastq.xz)
raw2=$(ls ../Raw/190822_AHFN3KDRXX-Lane1/*_S${SLURM_ARRAY_TASK_ID}_*R2*.fastq.xz)

# Cutadapt output
trimmed1=$(ls 1-trim/*_S${SLURM_ARRAY_TASK_ID}_*R1*.fastq.xz)
trimmed2=$(ls 1-trim/*_S${SLURM_ARRAY_TASK_ID}_*R2*.fastq.xz)

# As of FastQC 0.11.9 there is no support for xz yet, so use pipes
# FastQC requires two input filenames, so use named pipes
raw_fifo1=/tmp/$(basename ${raw1%.xz})
raw_fifo2=/tmp/$(basename ${raw2%.xz})
trimmed_fifo1=/tmp/$(basename ${trimmed1%.xz})
trimmed_fifo2=/tmp/$(basename ${trimmed2%.xz})

rm -f $raw_fifo1 $raw_fifo2 $trimmed_fifo1 $trimmed_fifo2
mkfifo $raw_fifo1 $raw_fifo2 $trimmed_fifo1 $trimmed_fifo2

xzcat $raw1 > $raw_fifo1 &
xzcat $raw2 > $raw_fifo2 &
xzcat $trimmed1 > $trimmed_fifo1 &
xzcat $trimmed2 > $trimmed_fifo2 &

# Background the first so we can use 4 cores per job
fastqc --threads=2 $raw_fifo1 $raw_fifo2 -o 2-qc/Raw &
fastqc --threads=2 $trimmed_fifo1 $trimmed_fifo2 -o 2-qc/Trimmed

rm -f $raw_fifo1 $raw_fifo2 $trimmed_fifo1 $trimmed_fifo2
