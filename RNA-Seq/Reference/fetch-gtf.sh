#!/bin/sh -e

fetch="$(./find-fetch.sh)"

# GTF
# Can't guarantee this file or the chromosome files will always be available.
# You may need to edit this.
release=$(./reference-release)
gtf=Mus_musculus.GRCm38.$release.gtf
if [ ! -e $gtf ]; then
    $fetch ftp://ftp.ensembl.org/pub/release-$release/gtf/mus_musculus/$gtf.gz
    gunzip $gtf.gz
else
    printf "$gtf already exists.\nRemove it and run $0 again to replace.\n"
fi

