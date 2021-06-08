#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on raw and trimmed data for comparison
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#   History:
#   Date        Name        Modification
#   2019-09-13  Jason Bacon Begin
##########################################################################

##########################################################################
#   Main
##########################################################################

# SLURM parameters
#SBATCH --array=1-9 --cpus-per-task=4
#SBATCH --mem=1g
#SBATCH --output=Data/2-qc/slurm-%A_%a.out
#SBATCH --error=Data/2-qc/slurm-%A_%a.err

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
    uname -a > Data/2-qc/os-version.txt 2>&1
    fastqc --version > Data/2-qc/fastqc-version.txt 2>&1
fi

mkdir -p Data/2-qc/Raw Data/2-qc/Trimmed

# One iteration if running under SLURM, all iterations otherwise
printf "Running task $SLURM_ARRAY_TASK_ID...\n"

# Raw files from sequencing center
raw1=$(ls ../Raw/190822_AHFN3KDRXX-Lane2/*_S${SLURM_ARRAY_TASK_ID}_*R1*.fastq.xz)
raw2=$(ls ../Raw/190822_AHFN3KDRXX-Lane2/*_S${SLURM_ARRAY_TASK_ID}_*R2*.fastq.xz)

# Cutadapt output
trimmed1=$(ls Data/1-trim/*_S${SLURM_ARRAY_TASK_ID}_*R1*.fastq.xz)
trimmed2=$(ls Data/1-trim/*_S${SLURM_ARRAY_TASK_ID}_*R2*.fastq.xz)

base_raw1=$(basename ${raw1%.fastq.xz})
base_raw2=$(basename ${raw2%.fastq.xz})
base_trimmed1=$(basename ${trimmed1%.fastq.xz})
base_trimmed2=$(basename ${trimmed2%.fastq.xz})

# Background the first three to run 4 jobs in parallel
xzcat $raw1 | fastqc --threads=2 -o Data/2-qc/Raw stdin:$base_raw1 &
xzcat $raw2 | fastqc --threads=2 -o Data/2-qc/Raw stdin:$base_raw2 &
xzcat $trimmed1 | fastqc --threads=2 -o Data/2-qc/Raw stdin:$base_trimmed1 &
xzcat $trimmed2 | fastqc --threads=2 -o Data/2-qc/Raw stdin:$base_trimmed2
