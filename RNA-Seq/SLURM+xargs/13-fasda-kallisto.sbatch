#!/bin/sh -e

#
#   Dependencies:
#       Requires kallisto abundances.  Run after *-kallisto-quant.sbatch.

#SBATCH --output=Logs/13-fasda-kallisto/slurm-%A_%a.out
#SBATCH --error=Logs/13-fasda-kallisto/slurm-%A_%a.err

# Document software versions used for publication
uname -a
fasda --version
pwd

kallisto_dir=Results/09-kallisto-quant
fasda_dir=Results/13-fasda-kallisto

for cell_type in neuro chondro; do
    for condition in time1 time2 time3; do
	printf "Normalizing $condition...\n"
	time fasda normalize \
	    --output $fasda_dir/$cell_type-$condition-all-norm.tsv \
	    $kallisto_dir/$cell_type-*-rep*-$condition/abundance.tsv
    done

    printf "Computing fold-change...\n"
    time fasda fold-change \
	--output $fasda_dir/$cell_type-time1-time2-FC.txt \
	$fasda_dir/$cell_type-time1-all-norm.tsv \
	$fasda_dir/$cell_type-time2-all-norm.tsv

    time fasda fold-change \
	--output $fasda_dir/$cell_type-time1-time3-FC.txt \
	$fasda_dir/$cell_type-time1-all-norm.tsv \
	$fasda_dir/$cell_type-time3-all-norm.tsv

    time fasda fold-change \
	--output $fasda_dir/$cell_type-time2-time3-FC.txt \
	$fasda_dir/$cell_type-time2-all-norm.tsv \
	$fasda_dir/$cell_type-time3-all-norm.tsv
done
