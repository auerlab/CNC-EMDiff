#!/bin/sh -e

build=$(../Common/genome-build.sh)
release=$(../Common/genome-release.sh)
fetch=$(../Common/find-fetch.sh)
gff=$(Reference/gff-filename.sh)

Reference/fetch-gff.sh
cd Data/3-reference
cat $gff | head -20 | tail -19 | awk '{ printf("%s\t%s\n", $2, $4); }' \
    | sort -n > chromosome-sizes.tsv
