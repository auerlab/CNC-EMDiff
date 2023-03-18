#!/bin/sh -e

##########################################################################
#   GFF is used by downstream analysis, such as peak classification
##########################################################################

fetch=$(../Common/find-fetch.sh)
release=$(../Common/genome-release.sh)
gff=$(Reference/gff-filename.sh)

# macOS zcat looks for .Z extension, while Linux does not have gzcat
zcat='gunzip -c'

##########################################################################
# Ensembl combined GFFs are sorted lexically by chromosome, while BAMs are
# sorted numerically.  Build our own GFF by concatenating individual
# chromosome GFFs in numeric order.  Resorting a GFF is complicated due
# to the hierarchical sort order (all gene components directly under the
# gene, etc).
##########################################################################

cd Results/07-reference
rm -f $gff
site=http://ftp.ensembl.org/pub/release-$release/gff3/mus_musculus

# Keep header from first GFF
file=Mus_musculus.GRCm38.98.chromosome.1.gff3.gz
if [ ! -e $file ]; then
    printf "Fetching $file...\n"
    $fetch $site/$file
fi
$zcat $file | egrep '^##gff|^#!' > $gff

# Concatenate the rest without the header
for chrom in $(seq 2 19); do
    file=Mus_musculus.GRCm38.98.chromosome.$chrom.gff3.gz
    if [ ! -e $file ]; then
	printf "Fetching $file...\n"
	$fetch $site/$file
    fi
    $zcat $file | egrep -v '^##[a-z]|^#!' >> $gff
done

