#!/bin/sh -e

##########################################################################
#   Script description:
#       Convert MACS2 peaks to fixed-size peaklets
#
#   Usage:
#       SLURM Cluster:
#           sbatch 9-process-peaks.sbatch
#       No cluster:
#           ./slurm-sim 9-process-peaks.sbatch 1 2 |& tee 1.log
#
#   History:
#       Based on the work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#   Date        Name        Modification
#   2020-02-02  Jason Bacon Begin
##########################################################################

##########################################################################
#   Main
##########################################################################

# SLURM parameters
#SBATCH --array=1-2
#SBATCH --mem=100
#SBATCH --output=9-process-peaks/slurm-%A_%a.out
#SBATCH --error=9-process-peaks/slurm-%A_%a.err

if [ -z $SLURM_ARRAY_TASK_ID ]; then
    printf "$0 must be run under sbatch or ./slurm-sim\n"
    exit 1
fi

##########################################################################
#   In MACS2 peak call output, filter for high-confidence peaks,
#   generate peaklets around summits, and merge overlapping peaklets.
##########################################################################

# Document software versions used for publication
uname -a > 9-process-peaks/os-version.txt 2>&1
bedtools --version > 9-process-peaks/bedtools-version.txt 2>&1

# One iteration if running under SLURM, all iterations otherwise
printf "Processing task $task_id...\n"
case $task_id in
1)
    cell_type=CCA
    ;;

2)
    cell_type=NCA
    ;;

*)
    printf "Error: task_id must be 1 or 2.\n"
    exit 1

esac

cd 9-process-peaks
macs_dir=../7-macs-peaklets

##########################################################################
#   Filter for high-confidence peaks: p-value < 1e-10
##########################################################################

# https://github.com/taoliu/MACS
# Col 8 or narrow peaks file is -log_10 p-value.
# We want p-values <= 10^-10, so col 8 >= 10
p_value=10

awk -v p_value=$p_value '$8 >= p_value { print $0 }' \
    $macs_dir/ATAC-$cell_type-macs.peaklets_peaks.narrowPeak \
    > high-confidence-p$p_value-$cell_type-narrow.bed

# Remove same lines from summits file that were filtered out of narrowPeak
# FIXME: This can all be done with just the narrow peaks file
#        The summits file is redundant
awk -f ../filter.awk \
    -v file_to_filter=$macs_dir/ATAC-$cell_type-macs.peaklets_summits.bed \
    high-confidence-p$p_value-$cell_type-narrow.bed \
    > high-confidence-p$p_value-$cell_type-summits.bed

##########################################################################
#   Create intervals from summit +/-250 nt
##########################################################################

# https://bedtools.readthedocs.io/en/latest/content/tools/slop.html
# explains nicely why you should use bedtools rather than just awk this.
# 
#awk '{ printf("%s\t%s\t%s\t%s\t%s\n", $1, $2-250, $2+250, $4, $5) }' \
#    ATAC-$cell_type-macs.peaklets_summits.bed \
#    > ATAC-$cell_type-macs.peaklets-501.bed

bedtools slop -b 250 -i high-confidence-p$p_value-$cell_type-summits.bed \
    -g ../../RNA-Seq/Reference/chromosome-sizes.tsv \
    > high-confidence-p$p_value-$cell_type-501.bed

head -4 high-confidence-p$p_value-$cell_type-summits.bed \
	high-confidence-p$p_value-$cell_type-501.bed

##########################################################################
#   Merge overlapping peaklets
##########################################################################
S1__S1_L01_L001-nodup-uniq.bam      Replicate 1, time 1
# https://bedtools.readthedocs.io/en/latest/content/tools/merge.html
merged=high-confidence-p$p_value-$cell_type-501-merged.bed
bedtools merge -i high-confidence-p$p_value-$cell_type-501.bed \
    > $merged
head -4 $merged
# Add 3 dummy fields ". . *" to match format of Andrea's merged BED files
awk '{ printf("%s\t%s\t.\t*\n", $0, "chr"$1"_"$2"_"$3 ) }' $merged > $merged.tmp
mv -f $merged.tmp $merged

trintf "\nWord counts:\n"
wc $macs_dir/ATAC-$cell_type-macs.peaklets_summits.bed \
    high-confidence-p$p_value-$cell_type-summits.bed \
    high-confidence-p$p_value-$cell_type-501-merged.bed

rm -f high-confidence-p$p_value-$cell_type-narrow.bed \
    high-confidence-p$p_value-$cell_type-summits.bed \
    high-confidence-p$p_value-$cell_type-501.bed
