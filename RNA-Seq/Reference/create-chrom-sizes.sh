#!/bin/sh -e

build=$(../Common/genome-build.sh)
release=$(../Common/genome-release.sh)
fetch=$(../Common/find-fetch.sh)

cd Data/3-reference
gff=Mus_musculus.GRCm$build.$release.chr.gff3.gz
if [ ! -e $gff ]; then
    $fetch ftp://ftp.ensembl.org/pub/release-$release/gff3/mus_musculus/$gff
fi

gzcat $gff | head -20 | tail -19 | awk '{ printf("%s\t%s\n", $2, $4); }' \
    | sort -n > chromosome-sizes.tsv
