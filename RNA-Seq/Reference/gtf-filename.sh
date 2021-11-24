#!/bin/sh -e

##########################################################################
#   This script should no longer be needed, since kallisto -gtf appears
#   to work with GFF3.  Use the equivalent gff3 script instead.
##########################################################################

build=$(../Common/genome-build.sh)
release=$(../Common/genome-release.sh)
echo Mus_musculus.GRCm$build.$release.chr.gtf
