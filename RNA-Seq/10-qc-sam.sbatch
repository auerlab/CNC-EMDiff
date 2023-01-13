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

# Set job array to number of samples.
#SBATCH --array=1-18
# FIXME: Check actual memory use
#SBATCH --mem=1g
#SBATCH --output=Logs/10-qc-sam/slurm-%A_%a.out
#SBATCH --error=Logs/10-qc-sam/slurm-%A_%a.err

: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
fastqc --version
pwd

# One iteration if running under SLURM, all iterations otherwise
printf "Running task $SLURM_ARRAY_TASK_ID...\n"

infile=Data/09-kallisto-quant/*sample$SLURM_ARRAY_TASK_ID*/pseudoalignments.bam
mkdir Data/10-qc-sam/$SLURM_ARRAY_TASK_ID
fastqc $infile -o Data/10-qc-sam/$SLURM_ARRAY_TASK_ID
