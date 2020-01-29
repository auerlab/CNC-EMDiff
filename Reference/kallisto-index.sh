#!/bin/sh -e

fetch="$(./find-fetch.sh)"
kallisto_reference=mus_musculus.tar.gz
if [ ! -e $kallisto_reference ]; then
    $fetch https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/$kallisto_reference
fi
tar ztf $kallisto_reference
