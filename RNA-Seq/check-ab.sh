#!/bin/sh -e

##########################################################################
#   Use this to understand why fasda and sleuth disagree about
#   significance of a FC
##########################################################################

usage()
{
    printf "Usage: $0 ID\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
printf "Time1:\n"
awk -v id=$1 '$1 == id { print $4 }' \
    Data/09-kallisto-quant/chondro-sample*-time1/abun*.tsv
printf "Time2:\n"
awk -v id=$1 '$1 == id { print $4 }' \
    Data/09-kallisto-quant/chondro-sample*-time2/abun*.tsv
