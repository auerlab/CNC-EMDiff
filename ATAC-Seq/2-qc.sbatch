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

# Temporary until all nodes are upgraded
#SBATCH --exclude=compute-[001-008,012]

#SBATCH --array=1-18
#SBATCH --mem=1g
#SBATCH --output=2-qc/slurm-%A_%a.out
#SBATCH --error=2-qc/slurm-%A_%a.err

# CentOS pkgsrc
if [ -e /sharedapps/pkg-2019Q2/bin/fastqc ]; then
    # sh is bash on CentOS, so this is OK
    source /etc/bashrc
    module load /sharedapps/pkg-2019Q2/etc/modulefiles/pkgsrc/2019Q2
fi

raw1="../Raw/190822_AHFN3KDRXX-Lane1/*_S${SLURM_ARRAY_TASK_ID}_*R1*.fastq.gz"
raw2="../Raw/190822_AHFN3KDRXX-Lane1/*_S${SLURM_ARRAY_TASK_ID}_*R2*.fastq.gz"

# Cutadapt output
trimmed1="1-trim/*_S${SLURM_ARRAY_TASK_ID}_*R1*.fastq.gz"
trimmed2="1-trim/*_S${SLURM_ARRAY_TASK_ID}_*R2*.fastq.gz"

printf "$raw1 $raw2 $trimmed1 $trimmed2\n"
mkdir -p 2-qc/Raw 2-qc/Trimmed

fastqc $raw1 -o 2-qc/Raw
fastqc $raw2 -o 2-qc/Raw
fastqc $trimmed1 -o 2-qc/Trimmed
fastqc $trimmed2 -o 2-qc/Trimmed
