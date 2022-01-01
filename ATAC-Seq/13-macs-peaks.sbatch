#!/bin/sh -e

##########################################################################
#   Script description:
#       Call peaks in alignments PE ATAC-seq aligned (BWA) reads
#
#   Usage:
#       SLURM Cluster:
#           sbatch 8-macs-peaks.sbatch
#       No cluster:
#           ../slurm-sim 8-macs-peaks.sbatch |& tee 1.log
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

# SLURM parameters
#SBATCH --array=1-2
#SBATCH --mem=4g
#SBATCH --output=8-macs-peaks/slurm-%A_%a.out
#SBATCH --error=8-macs-peaks/slurm-%A_%a.err

if [ -z $SLURM_ARRAY_TASK_ID ]; then
    printf "$0 must be run under sbatch or ../slurm-sim\n"
    exit 1
fi

# Document software versions used for publication
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
    uname -a > 8-macs-peaks/os-version.txt 2>&1
    macs2 --version > 8-macs-peaks/macs2-version.txt 2>&1
fi

case $SLURM_ARRAY_TASK_ID in
1)
    cell_type=CCA
    ;;

2)
    cell_type=NCA
    ;;

*)
    printf "Error: SLURM_ARRAY_TASK_ID must be 1 or 2.\n"
    exit 1

esac

printf "$cell_type: Peak calling, merging replicates across all time points, calling peaks...\n"
# "mm" is the compiled-in Mus Musculus genome size.  See MACS2 docs.
# BAMPE = BAM Paired End
# Would it make any difference to use the merged BAMs instead?
macs2 callpeak --nomodel \
    -t 6-remove-duplicates/$cell_type*_L001-nodup-uniq.bam \
    -f BAMPE -g mm --call-summits -n ATAC-$cell_type \
    --keep-dup all --outdir 8-macs-peaks
