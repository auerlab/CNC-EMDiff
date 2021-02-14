#!/bin/sh -e

##########################################################################
#   Script description:
#       Extract MAPQ scores from filtered BAMs and compute some stats
#       
#   History:
#   Date        Name        Modification
#   2021-02-14  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 file [file ...] | tee mapq-stats.tsv\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# -lt 1 ]; then
    usage
fi

for file in $@; do
    printf "$file...\n"
    samtools view $file | mawk '{ print $5 }' > $file.mapq
    alignments=$(cat $file.mapq | wc -l)
    printf "Alignments:\t%10u\n" $alignments
    for mapq in 1 5 10 20; do
	below=$(mawk -v mapq=$mapq '$1 < mapq' $file.mapq | wc -l)
	percent=$(printf "$below / $alignments * 100\nquit\n" | bc -l)
	p_incorrect=$(printf "e($mapq / -10 * l(10))\nquit\n" | bc -l)
	printf "MAPQ < $mapq:\t%10u\t%7.2f%%\tP(incorrect): %0.2f\n" \
	    $below $percent $p_incorrect
    done
done
