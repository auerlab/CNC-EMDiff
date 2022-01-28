#!/bin/sh -e

if [ $0 != "Reference/build-genome.sh" ]; then
    cat << EOM

$0 must be run as Reference/build-genome.sh.

EOM
    exit 1
fi

fetch=$(../Common/find-fetch.sh)
build=$(../Common/genome-build.sh)
release=$(../Common/genome-release.sh)
genome=$(Reference/genome-filename.sh)

# Chromosome files
mkdir -p Data/07-reference
cd Data/07-reference
for chromosome in $(seq 1 19); do
    file=Mus_musculus.GRCm$build.dna.chromosome.$chromosome.fa.gz
    if [ ! -e $file ]; then
	$fetch http://ftp.ensembl.org/pub/release-$release/fasta/mus_musculus/dna/$file
    fi
    chromosome=$((chromosome + 1))
done

if [ ! -e $genome ]; then
    printf "Concatenating chromosome FASTAs...\n"
    for chrom in $(seq 1 19); do
	printf "$chrom "
	zcat Mus_musculus.GRCm$build.dna.chromosome.$chrom.fa.gz >> $genome
    done
    printf "\n"
else
    printf "Using existing $genome...\n"
fi

if [ ! -e $genome.fai ]; then
    printf "Creating index $genome.fai...\n"
    samtools faidx $genome      # Speed up gffread
fi
