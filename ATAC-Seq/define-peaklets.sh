#!/bin/sh -e

##########################################################################
#   Script description:
#       In MACS2 peak call output, filter for high-confidence peaks,
#       generate peaklets around summits, and merge overlapping peaklets.
#
#   History:
#   Date        Name        Modification
#   2020-02-20  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 0 ]; then
    usage
fi

# https://bedtools.readthedocs.io/en/latest/content/tools/slop.html
# explains nicely why you should use bedtools rather than just awk this.
# 
#awk '{ printf("%s\t%s\t%s\t%s\t%s\n", $1, $2-250, $2+250, $4, $5) }' \
#    ATAC.nodup.unique.macs.peaklets_summits.bed \
#    > ATAC.nodup.unique.macs.peaklets-500.bed

cd 7-macs-peaklets

##########################################################################
#   Filter for high-confidence peaks: p-value < 1e-10
##########################################################################

awk '$8 >= 10 { print $0 }' ATAC.nodup.unique.macs.peaklets_peaks.narrowPeak \
    > high-confidence-narrow.bed
awk -f ../filter.awk \
    -v file_to_filter=ATAC.nodup.unique.macs.peaklets_summits.bed \
    high-confidence-narrow.bed > high-confidence-summits.bed

##########################################################################
#   Create intervals from summit +/-250 nt
##########################################################################

# https://bedtools.readthedocs.io/en/latest/content/tools/slop.html
bedtools slop -b 250 -i high-confidence-summits.bed \
    -g ../../RNA-Seq/Reference/chromosome-sizes.tsv \
    > high-confidence-500.bed
head -4 high-confidence-summits.bed high-confidence-500.bed

##########################################################################
#   Merge overlapping peaklets
##########################################################################

# https://bedtools.readthedocs.io/en/latest/content/tools/merge.html
bedtools merge -i high-confidence-500.bed > high-confidence-500-merged.bed
head -4 high-confidence-500-merged.bed

printf "\nWord counts:\n"
wc ATAC.nodup.unique.macs.peaklets_summits.bed high-confidence-summits.bed \
    high-confidence-500-merged.bed

read -p "Remove intermediate files? [y]/n " cleanup
if [ 0$cleanup != 0n ]; then
    rm -f high-confidence-narrow.bed high-confidence-summits.bed \
	high-confidence-500.bed
fi
