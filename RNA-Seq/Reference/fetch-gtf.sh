#!/bin/sh -e

fetch="$(./find-fetch.sh)"

# GTF
# Can't guarantee this file or the chromosome files will always be available.
# You may need to edit this.
gtf=Mus_musculus.GRCm38.98.gtf
if [ ! -e $gtf ]; then
    $fetch ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/$gtf.gz
    gunzip $gtf.gz
fi

