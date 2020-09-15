#!/bin/sh -e

files=""
chr=1
while [ $chr -le 19 ]; do
    files="$files Mus_musculus.GRCm38.dna.chromosome.$chr.fa.gz"
    chr=$(( $chr + 1 ))
done
autosome_file=Mus_musculus.GRCm38.dna.autosomes.fa
if [ ! -e $autosome_file ]; then
    printf "Concatenating $files to $autosome_file...\n"
    gzcat $files > $autosome_file
fi
