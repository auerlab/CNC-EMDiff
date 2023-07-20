#!/bin/sh -e


##########################################################################
#   Function description:
#       Pause until user presses return
##########################################################################

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}

# Raw counts
kallisto=Results/09-kallisto-quant/chondro-sample1-rep1-time1/abundance.tsv
hisat2=Results/19-fasda-abundance-hisat2/chondro-sample1-rep1-time1-abundance.tsv
for transcript in $(awk '{ print $1 }' $kallisto); do
    echo '==='
    printf "Kallisto: "
    grep $transcript $kallisto
    printf "Hisat2:   "
    grep $transcript $hisat2
done | more

# Normalized counts
kallisto=Results/13-fasda-kallisto/chondro-time1-all-norm.tsv
hisat2=Results/20-fasda-fc-hisat2/chondro-time1-all-norm.tsv
for transcript in $(awk '{ print $1 }' $kallisto); do
    echo '==='
    printf "Kallisto: "
    grep $transcript $kallisto
    printf "Hisat2:   "
    grep $transcript $hisat2
done | more

# Fold-changes
kallisto=Results/13-fasda-kallisto/chondro-time1-time2-FC.txt
hisat2=Results/20-fasda-fc-hisat2/chondro-time1-time2-FC.txt
for transcript in $(awk '{ print $1 }' $kallisto); do
    printf "\n"
    head -n 1 $kallisto
    grep $transcript $kallisto
    grep $transcript $hisat2
done | more
