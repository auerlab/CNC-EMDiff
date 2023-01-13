#!/bin/sh -e

##########################################################################
#   This is so fast, it's not even worth parallelizing
##########################################################################

#SBATCH --output=Logs/19-fasda-hisat2/slurm-%A_%a.out
#SBATCH --error=Logs/19-fasda-hisat2/slurm-%A_%a.err

# Document software versions used for publication
uname -a
fasda --version
pwd

hisat2_dir=Data/17-hisat2-align
fasda_dir=Data/19-fasda-hisat2

# Need to fetch gff3 for computing abundances
Reference/fetch-gff.sh
gff_filename=$(Reference/gff-filename.sh)
echo $gff_filename
reference_dir=Data/07-reference

##########################################################################
#   Compute abundances
##########################################################################

for cell_type in chondro neuro; do
    for time in $(seq 1 3); do
	for rep in $(seq 1 3); do
	    file=`ls $hisat2_dir/$cell_type-*rep$rep-time$time.bam`
	    ab=${file%.bam}-abundance.tsv
	    if [ ! -e $ab ]; then
		printf "Computing abundances for $condition replicate $r...\n"
		fasda abundance \
		    $reference_dir/$gff_filename $file
		head $ab
	    fi
	done
    done
done
exit

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
