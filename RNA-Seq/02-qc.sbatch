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
# and reverse raw and forward and reverse trimmed
#SBATCH --array=1-18
#SBATCH --cpus-per-task=4
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem=1g
#SBATCH --output=Logs/02-qc/slurm-%A_%a.out
#SBATCH --error=Logs/02-qc/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}

# Document software versions used for publication
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
    uname -a > Logs/02-qc/os-version-$SLURM_JOB_ID.txt 2>&1
    fastqc --version > Logs/02-qc/fastqc-version-$SLURM_JOB_ID.txt 2>&1
fi
hostname

mkdir -p Data/02-qc/Raw Data/02-qc/Trimmed

# Raw files from sequencing center
raw1=$(ls Data/Raw-renamed/*-sample${SLURM_ARRAY_TASK_ID}-*R1*.fastq.xz)
raw2=$(ls Data/Raw-renamed/*-sample${SLURM_ARRAY_TASK_ID}-*R2*.fastq.xz)

# Cutadapt output
trimmed1=$(ls Data/01-trim/*-sample${SLURM_ARRAY_TASK_ID}-*R1*.fastq.xz)
trimmed2=$(ls Data/01-trim/*-sample${SLURM_ARRAY_TASK_ID}-*R2*.fastq.xz)

# Filename stems for fastqc output
stem_raw1=$(basename ${raw1%.fastq.xz})
stem_raw2=$(basename ${raw2%.fastq.xz})
stem_trimmed1=$(basename ${trimmed1%.fastq.xz})
stem_trimmed2=$(basename ${trimmed2%.fastq.xz})

# Background the first three to run 4 jobs in parallel
xzcat $raw1 | fastqc -o Data/02-qc/Raw stdin:$stem_raw1 &
xzcat $raw2 | fastqc -o Data/02-qc/Raw stdin:$stem_raw2 &
xzcat $trimmed1 | fastqc -o Data/02-qc/Trimmed stdin:$stem_trimmed1 &
xzcat $trimmed2 | fastqc -o Data/02-qc/Trimmed stdin:$stem_trimmed2
