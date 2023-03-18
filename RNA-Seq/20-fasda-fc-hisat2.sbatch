#!/bin/sh -e

#
#   Dependencies:
#       Requires hisat2 abundances.
#       Run after *-fasda-abundance-hisat2.sbatch

#SBATCH --output=Logs/20-fasda-fc-hisat2/slurm-%A_%a.out
#SBATCH --error=Logs/20-fasda-fc-hisat2/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
fasda --version
samtools --version
pwd

hisat2_dir=Results/17-hisat2-align
fasda_dir=Results/20-fasda-fc-hisat2

# Need to fetch gff3 for computing abundances
gff_filename=$(Reference/gff-filename.sh)
reference_dir=Results/07-reference

##########################################################################
#   Compute abundances
##########################################################################

# Example input file: neuro-sample10-rep1-time1-abundance.tsv

for cell_type in neuro chondro; do
    for condition in time1 time2 time3; do
	printf "Normalizing $condition...\n"
	ls $hisat2_dir/$cell_type-*-rep*-$condition-abundance.tsv
	time fasda normalize \
	    --output $fasda_dir/$cell_type-$condition-all-norm.tsv \
	    $hisat2_dir/$cell_type-*-rep*-$condition-abundance.tsv
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
