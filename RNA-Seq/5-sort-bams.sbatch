#!/bin/sh -e

#SBATCH --array=1-18 --mem=1g
#SBATCH --output=5-sort-bams/slurm-%A_%a.out
#SBATCH --error=5-sort-bams/slurm-%A_%a.err

mkdir -p 5-sort-bams/$SLURM_ARRAY_TASK_ID

mv 5-sort-bams/$SLURM_ARRAY_TASK_ID/pseudoalignments-sorted.bam \
    5-sort-bams/pseudoalignments-sorted-$SLURM_ARRAY_TASK_ID.bam
exit

srun samtools sort \
    -o 5-sort-bams/pseudoalignments-sorted-$SLURM_ARRAY_TASK_ID.bam \
    4-kallisto-quant/$SLURM_ARRAY_TASK_ID/pseudoalignments.bam
