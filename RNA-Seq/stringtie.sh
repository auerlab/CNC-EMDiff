#!/bin/sh

kallisto=Results/09-kallisto-quant/chondro-sample1-rep1-time1/abundance.tsv
hisat2=Results/17-hisat2-align/chondro-sample1-rep1-time1.bam
gff=Results/07-reference/Mus_musculus.GRCm39.109.chr.gff3

if [ ! -e stringtie.out ]; then
    stringtie -e \
	-G $gff \
	$hisat2 \
	-A stringtie.out \
	-o stringtie.gtf
fi

# head 
more $kallisto
wc $kallisto stringtie.out
for transcript in $(cat stringtie.out | cut -f 1); do
    echo $transcript
    printf "Kallisto:\t"
    awk -v t=$transcript '$1 ~ t { print $4 }' $kallisto
    printf "Stringtie:\t"
    awk -v t=$transcript '$1 == t { print $7 }' stringtie.out
    echo '==='
done | more
