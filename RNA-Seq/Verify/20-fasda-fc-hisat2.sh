#!/bin/sh -e

# Must be run from parent parent of Verify dir for relative paths to work
if [ `basename $(pwd)` == Verify ]; then
    cd ..
fi

gff=$(Reference/gff-filename.sh)
echo $gff
awk '$3 ~ "RNA$|transcript$|gene_segment$"' Results/07-reference/$gff | wc -l
printf "\nCounts below should be 1 more than above.\n\n"
wc -l Results/20-fasda-fc-hisat2/*-FC.txt

# Add more sophisticated checks here if desired
