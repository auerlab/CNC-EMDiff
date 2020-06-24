#!/bin/sh -e

release=$(../../RNA-Seq/Reference/reference-release)
# Used 99 for first ATAC run
# release=99
site=ftp://ftp.ensembl.org/pub/release-$release/fasta/mus_musculus/dna
curl --continue-at - --remote-name $site/CHECKSUMS $site/README

chr=1
while [ $chr -le 19 ]; do
    file=Mus_musculus.GRCm38.dna.chromosome.$chr.fa.gz
    printf "Downloading $file...\n"
    curl --continue-at - --remote-name $site/$file
    chr=$(( $chr + 1 ))
done
