#!/bin/sh -e

#SBATCH --array=1-2
#SBATCH --mem=3g
#SBATCH --output=10-diff-anal/slurm-%A_%a.out
#SBATCH --error=10-diff-anal/slurm-%A_%a.err

# Spoof SLURM_ARRAY_TASK_ID if not set for testing outside SLURM env
if [ -z $SLURM_ARRAY_TASK_ID ]; then
    if [ $# = 1 ]; then
	export SLURM_ARRAY_TASK_ID=$1
    else
	cat << EOM

$0: You must run under SLURM or provide an array index
of 1 or 2 on the command-line.

EOM
	exit 1
    fi
fi

case $SLURM_ARRAY_TASK_ID in
1)
    cell_type=CCA
    ;;
2)
    cell_type=NCA
    ;;
esac

../10-diff-anal.R $cell_type
