#!/bin/sh -e

proper_name=Reference/cdna.sh
if [ $0 != "$proper_name" ]; then
    printf "$0 must be run as $proper_name\n"
    printf "from inside the Reference directory.\n"
    exit 1
fi

cd Data/3-reference

# Need GTF for kallisto quant --genomebam in any case
../../Reference/fetch-gtf.sh

fetch=$(../../../Common/find-fetch.sh)
build=$(../../../Common/genome-build.sh)
release=$(../../../Common/genome-release.sh)
awk=$(../../../Common/find-awk.sh)

# Can't guarantee this file will always be available.
# You may need to edit this.
cdna=Mus_musculus.GRCm$build.cdna.all.fa.gz
if [ ! -e $cdna ]; then
    $fetch ftp://ftp.ensembl.org/pub/release-$release/fasta/mus_musculus/cdna/$cdna
else
    printf "$cdna already exists.  Remove and rerun to replace.\n"
fi

# For compatibility with gtf_to_fasta alternative approach, don't gzip output
reference=$(../../Reference/reference-filename.sh)
set -x
zcat $cdna | $awk -f ../../Reference/remove-xy.awk > $reference
