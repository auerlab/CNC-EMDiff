#!/bin/sh -e

fetch=$(../../Common/find-fetch.sh)
release=$(../../Common/kallisto-release.sh)
kallisto_reference=mus_musculus.tar.gz
if [ ! -e $kallisto_reference ]; then
    $fetch https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-$release/$kallisto_reference
fi
tar ztf $kallisto_reference
