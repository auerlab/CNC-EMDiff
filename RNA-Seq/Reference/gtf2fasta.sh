#!/bin/sh -e

printf "FIXME: Update to run from Data/3-reference\n"
exit

proper_name="./gtf2fasta.sh"
if [ $0 != $proper_name ]; then
    printf "$0 must be run as $proper_name\n"
    printf "from inside the Reference directory.\n"
    exit 1
fi

./fetch-gtf.sh
fetch=$(../../Common/find-fetch.sh)
build=$(../../Common/genome-build.sh)
release=$(../../Common/genome-release.sh)

# Chromosome files
chromosome=1
while [ $chromosome -le 19 ]; do
    file=Mus_musculus.GRCm$build.dna.chromosome.$chromosome.fa.gz
    if [ ! -e $file ]; then
	$fetch http://ftp.ensembl.org/pub/release-$release/fasta/mus_musculus/dna/$file
    fi
    chromosome=$((chromosome + 1))
done

set -x
genome='all-but-xy.fa'
if [ ! -e $genome ]; then
    zcat Mus_musculus.GRCm$build.dna.chromosome.[0-9]*.gz > $genome
fi

# https://github.com/griffithlab/rnaseq_tutorial/wiki/Kallisto
transcripts=all-but-xy.transcripts.fa
if [ ! -e $transcripts ]; then
    gtf_to_fasta $gtf $genome $transcripts
fi

reference=$(./reference-filename.sh)
set -x
perl -ne 'if ($_ =~/^\>\d+\s+\w+\s+(ERCC\S+)[\+\-]/){print ">$1\n"}elsif($_ =~ /\d+\s+(ENST\d+)/){print ">$1\n"}else{print $_}' \
    $transcripts > $reference
