#!/bin/sh -e

fetch=$(../Common/find-fetch.sh)
release=$(../Common/genome-release.sh)
gff=$(Reference/gff-filename.sh)

# GFF
# Can't guarantee this file or the chromosome files will always be available.
# You may need to edit this.
cd Data/3-reference
if [ ! -e $gff.gz ]; then
    printf "Fetching $gff.gz...\n"
    $fetch ftp://ftp.ensembl.org/pub/release-$release/gff3/mus_musculus/$gff.gz
fi

if [ ! -e $gff ]; then
    # Filter for autosomes during decompress
    printf "Uncompressing and filtering $gff...\n"
    gunzip --stdout $gff.gz | awk '$1 ~ "^[0-9]"' > $gff
fi
