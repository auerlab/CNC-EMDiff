#!/bin/sh -e

proper_name="Reference/gff2fasta.sh"
if [ $0 != $proper_name ]; then
    printf "$0 must be run as $proper_name\n"
    printf "from inside the Reference directory.\n"
    exit 1
fi

Reference/fetch-gff.sh
fetch=$(../Common/find-fetch.sh)
build=$(../Common/genome-build.sh)
release=$(../Common/genome-release.sh)
reference=$(Reference/reference-filename.sh)
gff=$(Reference/gff-filename.sh)

cd Data/3-reference

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
    zcat Mus_musculus.GRCm$build.dna.chromosome.[0-9]*.gz > $genome
else
    printf "Using existing $genome...\n"
fi

# https://github.com/griffithlab/rnaseq_tutorial/wiki/Kallisto
transcripts=all-but-xy.transcripts.fa
if [ ! -e $transcripts ]; then
    # gtf_to_fasta is part of tophat, which is obsolete
    # gtf_to_fasta $gtf $genome $transcripts
    
    # Maybe?
    # bedtools getfasta -fi $genome -bed $gff > $transcripts
    
    # Recommended by Biostar RNA-Seq by Example
    # Warning: couldn't find fasta record for 'MT'!
    # Error: no genomic sequence available (check -g option!).
    # Abort trap (core dumped)
    printf "Converting $gff to $transcripts...\n"
    samtools faidx $genome      # Speed up gffread
    seq 1 19 > all-but-xy.list  # Mouse autosomes
    gffread --ids all-but-xy.list -w $transcripts -g $genome $gff
else
    printf "Using existing $transcripts...\n"
fi

# Tidy up headers
printf "Cleaning up headers, output in $reference...\n"
perl -ne 'if ($_ =~/^\>\d+\s+\w+\s+(ERCC\S+)[\+\-]/){print ">$1\n"}elsif($_ =~ /\d+\s+(ENST\d+)/){print ">$1\n"}else{print $_}' \
    $transcripts > $reference
