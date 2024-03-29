#!/bin/sh -e

##########################################################################
#   Script description:
#       Count unique differentially accessible peaks based on DESeq2
#       results.
#       
#   History:
#   Date        Name        Modification
#   2020-06-10  Jason Bacon Begin
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

# Record OS and software versions

# Use fgrep -v baseMean to strip headers if present
cd 10-diff-anal
for cell_type in CCA NCA; do

    printf "\n$cell_type LFC > 0:\n"
    rm -f peak-names.txt
    for time_interval in T1-vs-T0 T2-vs-T0 T2-vs-T1; do
	fgrep -v baseMean $cell_type-$time_interval.tsv \
	    | awk '$3 > 0 && $7 < 0.05 { print $1 }' >> peak-names.txt
    done
    printf "Sum across time intervals: "
    cat peak-names.txt | wc -l
    printf "Unique peaks:              "
    sort peak-names.txt | uniq | wc -l

    printf "\n$cell_type LFC < 0:\n"
    rm -f peak-names.txt
    for time_interval in T1-vs-T0 T2-vs-T0 T2-vs-T1; do
	fgrep -v baseMean $cell_type-$time_interval.tsv \
	    | awk '$3 < 0 && $7 < 0.05 { print $1 }' >> peak-names.txt
    done
    printf "Sum across time intervals: "
    cat peak-names.txt | wc -l
    printf "Unique peaks:              "
    sort peak-names.txt | uniq | wc -l
done
rm -f peak-names.txt
