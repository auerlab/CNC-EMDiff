#!/bin/sh -e

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}


# Must be run from parent parent of Verify dir for relative paths to work
if [ `basename $(pwd)` == Verify ]; then
    cd ..
fi

gff=$(Reference/gff-filename.sh)
echo $gff
awk '$3 ~ "RNA$|transcript$|gene_segment$"' Results/07-reference/$gff | wc -l
printf "\nCounts below should be 1 more than above.\n\n"
wc -l Results/19-fasda-abundance-hisat2/*.tsv
pause

kallisto=Results/09-kallisto-quant/chondro-sample1-rep1-time1/abundance.tsv
wc -l $kallisto
hisat2=Results/19-fasda-abundance-hisat2/chondro-sample1-rep1-time1-abundance.tsv
printf "%20s%20s%20s\n" "Feature ID" "Kallisto" "Hisat2"
for id in $(awk '{ print $1 }' $kallisto | fgrep -v target_id | head -40); do
    kc=$(awk -v id=$id '$1 ~ id { print $4 }' $kallisto)
    hc=$(awk -v id=$id '$1 ~ id { print $4 }' $hisat2)
    printf "%20s%20s%20s\n" $id $kc $hc
done | more

# Add more sophisticated checks here if desired
