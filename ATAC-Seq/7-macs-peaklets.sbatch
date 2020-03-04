#!/usr/bin/env bash

## Make peaks from PE ATAC-seq aligned (BWA) reads using MACS2
## Then we can quantify peaks using featureCounts

# Mem requirement observed on Peregrine
#SBATCH --array=1-2
#SBATCH --mem=4g
#SBATCH --output=7-macs-peaklets/slurm-%A_%a.out
#SBATCH --error=7-macs-peaklets/slurm-%A_%a.err

case $SLURM_ARRAY_TASK_ID in
1)
    cell_type=CCA
    ;;

2)
    cell_type=NCA
    ;;

*)
    printf "Error: SLURM_ARRAY_TASK_ID must be 1 or 2.\n"
    exit 1

esac

echo Peak calling, merging replicates across all time points, calling peaklets...
# "mm" is the compiled-in Mus Musculus genome size.  See MACS2 docs.
# BAMPE = BAM Paired End
# Does order of the filenames matter?
macs2 callpeak --nomodel -t 4-bwa-mem/$cell_type*_L001-nodup-uniq.bam \
    -f BAMPE -g mm --call-summits -n ATAC-$cell_type-macs.peaklets \
    --keep-dup all --outdir 7-macs-peaklets

##########################################################################
#   In MACS2 peak call output, filter for high-confidence peaks,
#   generate peaklets around summits, and merge overlapping peaklets.
##########################################################################

cd 7-macs-peaklets

##########################################################################
#   Filter for high-confidence peaks: p-value < 1e-10
##########################################################################

# https://github.com/taoliu/MACS
# Col 8 or narrow peaks file is -log_10 p-value.
# We want p-values <= 10^-10, so col 8 >= 10
awk '$8 >= 10 { print $0 }' ATAC-$cell_type-macs.peaklets_peaks.narrowPeak \
    > high-confidence-$cell_type-narrow.bed

# Remove same lines from summits file that were filtered out of narrowPeak
awk -f ../filter.awk \
    -v file_to_filter=ATAC-$cell_type-macs.peaklets_summits.bed \
    high-confidence-$cell_type-narrow.bed \
    > high-confidence-$cell_type-summits.bed

##########################################################################
#   Create intervals from summit +/-250 nt
##########################################################################

# https://bedtools.readthedocs.io/en/latest/content/tools/slop.html
# explains nicely why you should use bedtools rather than just awk this.
# 
#awk '{ printf("%s\t%s\t%s\t%s\t%s\n", $1, $2-250, $2+250, $4, $5) }' \
#    ATAC-$cell_type-macs.peaklets_summits.bed \
#    > ATAC-$cell_type-macs.peaklets-501.bed

bedtools slop -b 250 -i high-confidence-$cell_type-summits.bed \
    -g ../../RNA-Seq/Reference/chromosome-sizes.tsv \
    > high-confidence-$cell_type-501.bed

head -4 high-confidence-$cell_type-summits.bed \
	high-confidence-$cell_type-501.bed

##########################################################################
#   Merge overlapping peaklets
##########################################################################

# https://bedtools.readthedocs.io/en/latest/content/tools/merge.html
bedtools merge -i high-confidence-$cell_type-501.bed \
    > high-confidence-$cell_type-501-merged.bed
head -4 high-confidence-$cell_type-501-merged.bed

printf "\nWord counts:\n"
wc ATAC-$cell_type-macs.peaklets_summits.bed \
    high-confidence-$cell_type-summits.bed \
    high-confidence-$cell_type-501-merged.bed

rm -f high-confidence-$cell_type-narrow.bed \
    high-confidence-$cell_type-summits.bed \
    high-confidence-$cell_type-501.bed
