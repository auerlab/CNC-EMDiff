#!/bin/sh -e

site=ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna
curl --continue-at - --remote-name $site/CHECKSUMS $site/README

chr=1
while [ $chr -le 22 ]; do
    for variant in "" "_rm" "_sm"; do
	file=Mus_musculus.GRCm38.dna$variant.chromosome.$chr.fa.gz
	printf "Downloading $file...\n"
	curl --continue-at - --remote-name $site/$file
    done
    chr=$(( $chr + 1 ))
done
