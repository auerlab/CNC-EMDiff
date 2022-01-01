#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on raw and trimmed data for comparison
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#
#   History:
#   Date        Name        Modification
#   2019-09-13  Jason Bacon Begin
##########################################################################

# Each invocation of this script will run 4 fastqc processes, forward
# and reverse
#SBATCH --array=1-18
#SBATCH --cpus-per-task=2
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem=1g
#SBATCH --output=Logs/04-qc/slurm-%A_%a.out
#SBATCH --error=Logs/04-qc/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
fastqc --version
pwd

mkdir -p Data/04-qc

# Cutadapt output
trimmed1=$(ls Data/03-trim/*-sample${SLURM_ARRAY_TASK_ID}-*R1*.fastq.xz)
trimmed2=$(ls Data/03-trim/*-sample${SLURM_ARRAY_TASK_ID}-*R2*.fastq.xz)

# Filename stems for fastqc output
stem_trimmed1=$(basename ${trimmed1%.fastq.xz})
stem_trimmed2=$(basename ${trimmed2%.fastq.xz})

# Background the first three to run 4 jobs in parallel
xzcat $trimmed1 | fastqc -o Data/04-qc stdin:$stem_trimmed1 &
xzcat $trimmed2 | fastqc -o Data/04-qc stdin:$stem_trimmed2
