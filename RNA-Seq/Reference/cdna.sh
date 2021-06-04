#!/bin/sh -e

#SBATCH --ntasks=1 --mem=1g

proper_name=./cdna.sh
if [ $0 != "$proper_name" ]; then
    printf "$0 must be run as $proper_name\n"
    printf "from inside the Reference directory.\n"
    exit 1
fi

# Need GTF for kallisto quant --genomebam in any case
./fetch-gtf.sh

fetch=$(../../Common/find-fetch.sh)
build=$(../../Common/genome-build.sh)
release=$(../../Common/genome-release.sh)

# Can't guarantee this file will always be available.
# You may need to edit this.
cdna=Mus_musculus.GRCm$build.cdna.all.fa.gz
if [ ! -e $cdna ]; then
    $fetch ftp://ftp.ensembl.org/pub/release-$release/fasta/mus_musculus/cdna/$cdna
else
    printf "$cdna already exists.  Remove and rerun to replace.\n"
fi

# For compatibility with gtf_to_fasta alternative approach, don't gzip output
reference=$(./reference-filename.sh)
set -x
zcat $cdna | awk -f remove-xy.awk > $reference
