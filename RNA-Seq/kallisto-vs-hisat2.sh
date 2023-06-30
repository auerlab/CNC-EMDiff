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

# Abundances
kallisto=Results/13-fasda-kallisto/chondro-time1-all-norm.tsv
hisat2=Results/20-fasda-fc-hisat2/chondro-time1-all-norm.tsv
for transcript in $(awk '{ print $1 }' $kallisto | head -20); do
    echo '==='
    printf "Kallisto: "
    grep $transcript $kallisto
    printf "Hisat2:   "
    grep $transcript $hisat2
done
exit

# Fold-changes
kallisto=Results/13-fasda-kallisto/chondro-time1-time2-FC.txt
hisat2=Results/20-fasda-fc-hisat2/chondro-time1-time2-FC.txt
for transcript in $(awk '{ print $1 }' $kallisto | head -20); do
    echo $transcript
    grep $transcript $kallisto
    grep $transcript $hisat2
done
