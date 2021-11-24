#!/bin/sh -ex

build=$(../Common/genome-build.sh)
release=$(../Common/genome-release.sh)
fetch=$(../Common/find-fetch.sh)
gff=$(Reference/gff-filename.sh)

Reference/fetch-gff.sh
cd Data/3-reference
# GFF col 1 = chromosome name
# col 4 = start pos, always 1 for a chromosome feature
# col 5 = end pos = size of chromosome
awk '$3 == "chromosome" { printf("%s\t%s\n", $1, $5); }' $gff \
    | sort -n > chromosome-sizes.tsv
cat chromosome-sizes.tsv
