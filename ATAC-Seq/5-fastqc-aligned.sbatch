#!/bin/sh -e

#SBATCH --array=1-18
# FIXME: Check this
#SBATCH --mem=1g

# Fake SLURM_ARRAY_TASK_ID for testing script outside SLURM env
: ${SLURM_ARRAY_TASK_ID:=1}

cd 5-qc-sam
infile=../4-BWA-mem/*_S${SLURM_ARRAY_TASK_ID}_L001.sam
fastqc $infile
