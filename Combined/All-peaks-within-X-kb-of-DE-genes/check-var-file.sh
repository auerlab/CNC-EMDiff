#!/bin/sh -e

##########################################################################
#   Script description:
#       Show stats on variable-column output files to verify gene counts
#       
#   History:
#   Date        Name        Modification
#   2020-09-24  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 gene-peak-file\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi

file=$1
var_file=${file%.tsv}-var-cols.tsv

total_genes=0
for cols in $(seq 12 2 40); do
    lines=$(awk '{ print NF }' $var_file | fgrep $cols | wc -l)
    genes=$(( (cols - 10) / 2 ))
    printf "Peaks near %2u genes: $lines\n" $genes
    total_genes=$(( total_genes + lines * genes ))
done

multi_gene=$(awk 'NF > 12' $var_file | wc -l)
printf "\n"
printf "Total peaks near multiple genes:     %u\n" $multi_gene

printf "Total genes in variable-column file: %u\n" $total_genes

fixed_peaks=$(cat $file | wc -l)
printf "Total genes in fixed-column file:    %u\n" $fixed_peaks
