#!/bin/sh -e

##########################################################################
#   Script description:
#       QC kallisto alignments
#
#   Usage:
#       SLURM cluster:
#           sbatch 5-qc.sbatch
#       No cluster:
#           ../slurm-sim 5-qc.sbatch |& tee 5.log
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

# SLURM parameters
#SBATCH --array=1-18
# FIXME: Check actual memory use
#SBATCH --mem=1g
#SBATCH --output=5-qc-sam/slurm-%A_%a.out
#SBATCH --error=5-qc-sam/slurm-%A_%a.err

printf "FIXME: Adapt from ATAC-Seq script.\n"
exit

if [ -z $SLURM_ARRAY_TASK_ID ]; then
    printf "$0 must be run under sbatch or ../slurm-sim\n"
    exit 1
fi

# Document software versions used for publication
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
    uname -a > 5-qc-sam/os-version.txt 2>&1
    fastqc --version > 5-qc-sam/fastqc-version.txt 2>&1
fi

# One iteration if running under SLURM, all iterations otherwise
printf "Running task $SLURM_ARRAY_TASK_ID...\n"

infile=4-bwa-mem/*_S${SLURM_ARRAY_TASK_ID}_L001.sam
fastqc --version > 5-qc-sam/fastqc-version.txt 2>&1
fastqc $infile -o 5-qc-sam
