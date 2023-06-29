#!/bin/sh -e

build=$(../Common/genome-build.sh)
release=$(../Common/genome-release.sh)

echo Mus_musculus.GRCm$build.$release.chr.gff3
