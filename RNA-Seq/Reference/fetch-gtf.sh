#!/bin/sh -e

fetch=$(../../Common/find-fetch.sh)
gtf=$(./gtf-filename.sh)
release=$(../../Common/genome-release.sh)

# GTF
# Can't guarantee this file or the chromosome files will always be available.
# You may need to edit this.
if [ ! -e $gtf ]; then
    $fetch ftp://ftp.ensembl.org/pub/release-$release/gtf/mus_musculus/$gtf.gz
    gunzip $gtf.gz
else
    printf "$gtf already exists.\nRemove it and run $0 again to replace.\n"
fi

