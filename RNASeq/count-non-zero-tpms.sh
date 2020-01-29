#!/bin/sh -e

##########################################################################
#   Script description:
#       Count entries with non-zero TPM values in kallisto output
#       
#   History:
#   Date        Name        Modification
#   2019-10-09  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 kallisto-directory\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi

dir=$1
sample=1
while [ $sample -le 18 ]; do
    printf "$dir/$sample/abundance.tsv\t"
    awk '$5 != 0 { print $5 }' 4-kallisto-quant/$sample/abundance.tsv | wc -l
    sample=$((sample+1))
done
