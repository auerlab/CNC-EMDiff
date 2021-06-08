#!/bin/sh -e

##############################################################################
#   Find overlap of DA peaks from neuro and chondro that are within 250
#   kb of commonly regulated DE genes
#
#     1. report % overlapping peaks
#     2. report % overlapping peaks with matching temporal accessibility
#        patterns
#     3. for each overlap report DA temporal accessibility patterns for
#        chondro and neuro  2. Find overlap of DA peaks from neuro and chondro that are 
##############################################################################

# Filter peaks for differential accessibility across any time interval
# and genes in genes-common*.txt

##############################################################################
#
# Sample peak input:
#
# Chr     Start   End     Peak-name           
# 1       6272149 6272650 chr1-6272149-6272650
#
# T1-T0-LFC    T1-T0-APV   T2-T0-LFC    T2-T0-APV   T2-T1-LFC   T2-T1-APV
# -0.528       0.755       -0.227       0.689       0.301       0.533
#
# Gene-name Gene-Distance
# St18      214581
##############################################################################

peaks_dir=../All-peaks-within-X-kb-of-DE-genes/Proximal+distal-250kb
common_genes="$(cat Data/genes-common-*.txt | sort | uniq)"
common_genes=$(echo $common_genes)  # Convert newlines to spaces
echo $common_genes
printf "%u commonly regulated genes (union of up and down)\n" \
    $(echo $common_genes | wc -w)

for file in $peaks_dir/genes-CCA-c[1-3].bed-proximal+distal-peaks.tsv; do
    echo $file
    awk -v genestr="$common_genes" -f da+proximal-gene.awk $file | head
done
for file in $peaks_dir/genes-NCA-c[1-5].bed-proximal+distal-peaks.tsv; do
    echo $file
    awk -v genestr="$common_genes" -f da+proximal-gene.awk $file | head
done

