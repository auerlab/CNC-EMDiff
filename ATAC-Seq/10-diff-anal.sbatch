#!/bin/sh -e

#SBATCH --array=1-2
#SBATCH --mem=3g
#SBATCH --output=10-diff-anal/slurm-%A_%a.out
#SBATCH --error=10-diff-anal/slurm-%A_%a.err

# Spoof SLURM_ARRAY_TASK_ID if not set for testing outside SLURM env
: ${SLURM_ARRAY_TASK_ID:=1}

case $SLURM_ARRAY_TASK_ID in
1)
    cell_type=CCA
    ;;
2)
    cell_type=NCA
    ;;
esac

cd 10-diff-anal
../10-diff-anal.R $cell_type
