#!/bin/sh -e

##########################################################################
#   Script description:
#       Call peaks in alignments PE ATAC-seq aligned (BWA) reads
#
#   History:
#       Based on the work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#   Date        Name        Modification
#   2020-03-10  Jason Bacon Begin
##########################################################################

##########################################################################
#   Main
##########################################################################

usage()
{
    printf "Usage: $0 chondro|neuro\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
cell_type=$1

# Document software versions used for publication
uname -a
macs3 --version
pwd

printf "$cell_type: Peak calling, merging replicates across all time points, calling peaks...\n"
# "mm" is the compiled-in Mus Musculus genome size.  See MACS3 docs.
# BAMPE = BAM Paired End
# Would it make any difference to use the merged BAMs instead?
set -x
macs3 callpeak --nomodel \
    -t Results/10-remove-duplicates/$cell_type-*-mapq1.bam \
    -f BAMPE -g mm --call-summits -n ATAC-$cell_type \
    --keep-dup all --outdir Results/14-macs-peaks
