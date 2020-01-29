#!/bin/sh -e

##########################################################################
#   Script description:
#       Show all read lengths in sample files
#       
#   History:
#   Date        Name        Modification
#   2019-10-08  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 directory\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi

dir=$1
for file in $dir/*.gz; do
    echo $file
    zcat $file | awk 'NR % 4 == 2 { print length }' \
	| grep 20 #-v 101 | sort -nr | tee $file.read-lengths.txt
done
