#!/bin/sh -e

##########################################################################
#   Script description:
#       Convert MACS2 peaks to fixed-size peaklets
#
#   Usage:
#       SLURM Cluster:
#           sbatch 9-process-peaks.sbatch
#       No cluster:
#           ../slurm-sim 9-process-peaks.sbatch |& tee 1.log
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
    printf "$0 must be run under sbatch or ../slurm-sim\n"
    exit 1
fi

##########################################################################
#   In MACS2 peak call output, filter for high-confidence peaks,
#   generate peaklets around summits, and merge overlapping peaklets.
#   Centrino? requires that summits be in the middle of peaklets (Ava)
#   Some MEME tools require all peaklets to be the same size (Ava)
##########################################################################

# Document software versions used for publication
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
    uname -a > 9-process-peaks/os-version.txt 2>&1
    bedtools --version > 9-process-peaks/bedtools-version.txt 2>&1
fi

# One iteration if running under SLURM, all iterations otherwise
printf "Processing task $SLURM_ARRAY_TASK_ID...\n"
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

cd 9-process-peaks
macs_dir=../8-macs-peaks

##########################################################################
#   Filter for high-confidence peaks: p-value < 1e-10
##########################################################################

# https://github.com/taoliu/MACS
# Col 8 or narrow peaks file is -log_10 p-value.
# We want p-values <= 10^-10, so col 8 >= 10
p_value=10

awk -v p_value=$p_value '$8 >= p_value { print $0 }' \
    $macs_dir/ATAC-${cell_type}_peaks.narrowPeak \
    > p$p_value-$cell_type-narrow.bed

# Remove same lines from summits file that were filtered out of narrowPeak
# FIXME: This can all be done with just the narrow peaks file
#        The summits file is redundant
# Output format matches summits file: chr start end peak-id LFC
# LFC is always 0 here.  These peaks will be merged with DESeq2 output
# later to fill in LFC values for DA peaks.
awk -f ../Aux/filter.awk \
    -v file_to_filter=$macs_dir/ATAC-${cell_type}_summits.bed \
    p$p_value-$cell_type-narrow.bed \
    > p$p_value-$cell_type-summits.bed

##########################################################################
#   Create intervals from summit +/-250 nt
##########################################################################

# https://bedtools.readthedocs.io/en/latest/content/tools/slop.html
# explains nicely why you should use bedtools rather than just awk this.
# 
#awk '{ printf("%s\t%s\t%s\t%s\t%s\n", $1, $2-250, $2+250, $4, $5) }' \
#    ATAC-${cell_type}_summits.bed \
#    > ATAC-$cell_type-501.bed

bed_sort="sort -k1,1 -k2,2n"
bedtools slop -b 250 -i p$p_value-$cell_type-summits.bed \
    -g ../../RNA-Seq/Reference/chromosome-sizes.tsv \
    | $bed_sort > p$p_value-$cell_type-501.bed

head -4 p$p_value-$cell_type-summits.bed \
	p$p_value-$cell_type-501.bed

##########################################################################
#   Merge overlapping peaklets
#   This seems to prevent DiffBind and DESeq2 from messing with the
#   peaklet boundaries, so DESeq2 outputs will contain the same peaklets
#   as the file $merged here.
##########################################################################
# S1__S1_L01_L001-nodup-uniq.bam      Replicate 1, time 1
# https://bedtools.readthedocs.io/en/latest/content/tools/merge.html

merged=p$p_value-$cell_type-501-merged.bed
bedtools merge -i p$p_value-$cell_type-501.bed > $merged
head -4 $merged
# Add 3 dummy fields ". . *" to match format of Andrea's merged BED files
awk '{ printf("%s\t%s\t%s\t%s\t0\n", $1, $2, $3, "chr"$1"-"$2"-"$3 ) }' \
    $merged > $merged.tmp
mv -f $merged.tmp $merged

printf "\nWord counts:\n"
wc $macs_dir/ATAC-${cell_type}_summits.bed \
    p$p_value-$cell_type-summits.bed \
    p$p_value-$cell_type-501.bed

rm -f p$p_value-$cell_type-narrow.bed \
    p$p_value-$cell_type-summits.bed
