#!/bin/sh -e

fetch=$(../Common/find-fetch.sh)
release=$(../Common/genome-release.sh)
gff=$(Reference/gff-filename.sh)

# GFF
# Can't guarantee this file or the chromosome files will always be available.
# You may need to edit this.
cd Data/3-reference
if [ ! -e $gff ]; then
    $fetch ftp://ftp.ensembl.org/pub/release-$release/gff3/mus_musculus/$gff.gz
    gunzip $gff.gz
fi
