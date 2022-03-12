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

#SBATCH --array=1-2
#SBATCH --mem=4g
#SBATCH --output=Logs/13-macs-peaks/slurm-%A_%a.out
#SBATCH --error=Logs/13-macs-peaks/slurm-%A_%a.err

: ${SLURM_ARRAY_TASK_ID:=1}

uname -a
macs2 --version
pwd

case $SLURM_ARRAY_TASK_ID in
1)
    cell_type=chondro
    ;;

2)
    cell_type=neuro
    ;;

*)
    printf "Error: SLURM_ARRAY_TASK_ID must be 1 or 2.\n"
    exit 1

esac

printf "$cell_type: Peak calling, merging replicates across all time points, calling peaks...\n"
# "mm" is the compiled-in Mus Musculus genome size.  See MACS2 docs.
# BAMPE = BAM Paired End
# Would it make any difference to use the merged BAMs instead?
set -x
macs2 callpeak --nomodel \
    -t Data/09-remove-duplicates/$cell_type-*-mapq1.bam \
    -f BAMPE -g mm --call-summits -n ATAC-$cell_type \
    --keep-dup all --outdir Data/13-macs-peaks
