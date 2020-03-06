#!/bin/sh -e

##########################################################################
#   Prepare input file for loading into R for DESeq2
##########################################################################

# BAM
awk '{ printf("%s\t%s\t%s\t%s\n", $1, $2, $3, "chr"$1"-"$2"-"$3) }' \
    7-macs-peaklets/high-confidence-CCA-500-merged.bed
