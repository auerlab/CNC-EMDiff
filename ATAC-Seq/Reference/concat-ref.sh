#!/bin/sh -e

for variant in "" "_rm" "_sm"; do
    files=""
    chr=1
    while [ $chr -le 22 ]; do
	files="$files Mus_musculus.GRCm38.dna$variant.chromosome.$chr.fa.gz"
	chr=$(( $chr + 1 ))
    done
    autosome_file=Mus_musculus.GRCm28.dna$variant.autosomes.fa.gz
    printf "Concatenating $files to $autosome_file...\n"
    gzcat $files | gzip - > $autosome_file
done
