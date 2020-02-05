#!/bin/sh -e

#SBATCH --ntasks=1 --mem=1g

proper_name="./gtf2fasta.sh"
if [ $0 != $proper_name ]; then
    printf "$0 must be run as $proper_name\n"
    printf "from inside the Reference directory.\n"
    exit 1
fi

fetch="$(./find-fetch.sh)"

# GTF
# Can't guarantee this file or the chromosome files will always be available.
# You may need to edit this.
gtf=Mus_musculus.GRCm38.98.gtf
if [ ! -e $gtf ]; then
    $fetch ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/$gtf.gz
    gunzip $gtf.gz
fi

# Chromosome files
chromosome=1
while [ $chromosome -le 19 ]; do
    file=Mus_musculus.GRCm38.dna.chromosome.$chromosome.fa.gz
    if [ ! -e $file ]; then
	$fetch ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/$file
    fi
    chromosome=$((chromosome + 1))
done

set -x
genome='all-but-xy.fa'
if [ ! -e $genome ]; then
    zcat Mus_musculus.GRCm38.dna.chromosome.[0-9]*.gz > $genome
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
