#!/bin/sh -e

#SBATCH --array=1-18 --mem=1g

srun samtools sort \
    -o 4-kallisto-quant/$SLURM_ARRAY_TASK_ID/pseudoalignments-sorted.bam
    4-kallisto-quant/$SLURM_ARRAY_TASK_ID/pseudoalignments.bam
