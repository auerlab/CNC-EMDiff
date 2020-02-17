#!/bin/sh -e

#SBATCH --array=1-18
# FIXME: Check memory use
#SBATCH --mem=1g
#SBATCH --output=5-qc-sam/slurm-%A_%a.out
#SBATCH --error=5-qc-sam/slurm-%A_%a.err

# Spoof SLURM_ARRAY_TASK_ID if not set for testing outside SLURM env
: ${SLURM_ARRAY_TASK_ID:=1}

infile=4-bwa-mem/*_S${SLURM_ARRAY_TASK_ID}_L001.sam
fastqc $infile -o 5-qc-sam
