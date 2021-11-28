#!/bin/sh -e

proper_name="Reference/gtf2fasta.sh"
if [ $0 != $proper_name ]; then
    printf "$0 must be run as $proper_name\n"
    printf "from inside the Reference directory.\n"
    exit 1
fi

Reference/fetch-gtf.sh
fetch=$(../Common/find-fetch.sh)
build=$(../Common/genome-build.sh)
release=$(../Common/genome-release.sh)
transcriptome=$(Reference/transcriptome-filename.sh)
gtf=$(Reference/gtf-filename.sh)

cd Data/03-reference

# Chromosome files
chromosome=1
while [ $chromosome -le 19 ]; do
    file=Mus_musculus.GRCm$build.dna.chromosome.$chromosome.fa.gz
    if [ ! -e $file ]; then
	$fetch http://ftp.ensembl.org/pub/release-$release/fasta/mus_musculus/dna/$file
    fi
    chromosome=$((chromosome + 1))
done

genome='all-but-xy.fa'
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

# https://github.com/griffithlab/rnaseq_tutorial/wiki/Kallisto
if [ ! -e $transcriptome ]; then
    # gtf_to_fasta is part of tophat, which is obsolete
    # gtf_to_fasta $gtf $genome $transcriptome
    
    # Maybe?
    # bedtools getfasta -fi $genome -bed $gtf > $transcriptome
    
    # Recommended by Biostar RNA-Seq by Example
    # Warning: couldn't find fasta record for 'MT'!
    # Error: no genomic sequence available (check -g option!).
    # Abort trap (core dumped)
    printf "Converting $gtf to $transcriptome...\n"
    gffread -w $transcriptome -g $genome $gtf
else
    printf "Using existing $transcriptome...\n"
fi

# Tidy up headers: Actually has no effect
# perl -ne 'if ($_ =~/^\>\d+\s+\w+\s+(ERCC\S+)[\+\-]/){print ">$1\n"}elsif($_ =~ /\d+\s+(ENST\d+)/){print ">$1\n"}else{print $_}' \

