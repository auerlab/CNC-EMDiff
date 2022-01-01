#!/bin/sh -e

outfile=CCA-NCA-intersect-all.tsv
cca_in=9-process-peaks/p10-CCA-501-merged.bed
nca_in=9-process-peaks/p10-NCA-501-merged.bed
bedtools intersect \
    -a $cca_in \
    -b $nca_in \
    -wo | cut -f 4,9,11 | sort -nk 3 > $outfile

printf "CCA peaks:                 "
cat $cca_in | wc -l
printf "NCA peaks:                 "
cat $nca_in | wc -l
printf "Total overlaps:            "
cat $outfile | wc -l
for ov in 10 100 200 300 400 500 1000 10000; do
    printf "Overlaps < %5s bases:    " $ov
    awk -v ov=$ov '$3 < ov' $outfile | wc -l
done
