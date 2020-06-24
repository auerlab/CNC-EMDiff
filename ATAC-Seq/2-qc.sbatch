#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on raw and trimmed data for comparison
#
#   Usage:
#       SLURM cluster:
#           sbatch 2-qc.sbatch
#       No cluster:
#           ./slurm-sim 2-qc.sbatch 1 18 |& tee 2.log
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
#SBATCH --array=1-18
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

if [ -z $SLURM_ARRAY_TASK_ID ]; then
    printf "$0 must be run under sbatch or ./slurm-sim\n"
    exit 1
fi

# Defaults if not set by SLURM or user
: ${SLURM_CPUS_PER_TASK:=4}

# Document version of software used for publication
fastqc --version > 2-qc/fastqc-version.txt 2>&1

mkdir -p 2-qc/Raw 2-qc/Trimmed

# One iteration if running under SLURM, all iterations otherwise
printf "Running task $SLURM_ARRAY_TASK_ID...\n"

# UWM CentOS cluster
if [ -e /sharedapps/pkg-2019Q2/bin/fastqc ]; then
    # sh is bash on CentOS, so this is OK
    source /etc/bashrc
    module load /sharedapps/pkg-2019Q2/etc/modulefiles/pkgsrc/2019Q2
fi

# Raw files from sequencing center
raw1="../Raw/190822_AHFN3KDRXX-Lane1/*_S${SLURM_ARRAY_TASK_ID}_*R1*.fastq.gz"
raw2="../Raw/190822_AHFN3KDRXX-Lane1/*_S${SLURM_ARRAY_TASK_ID}_*R2*.fastq.gz"

# Cutadapt output
trimmed1="1-trim/*_S${SLURM_ARRAY_TASK_ID}_*R1*.fastq.gz"
trimmed2="1-trim/*_S${SLURM_ARRAY_TASK_ID}_*R2*.fastq.gz"

fastqc $raw1 -o 2-qc/Raw
fastqc $raw2 -o 2-qc/Raw
fastqc $trimmed1 -o 2-qc/Trimmed
fastqc $trimmed2 -o 2-qc/Trimmed
