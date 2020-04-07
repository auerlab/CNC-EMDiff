#!/bin/sh -e

#SBATCH --array=1-2
#SBATCH --mem=3g

case $SLURM_ARRAY_TASK_ID in
1)
    cell_type=CCA
    ;;
2)
    cell_type=NCA
    ;;
esac

Rscript 10-diff-anal.R $cell_type
