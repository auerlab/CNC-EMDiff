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
#SBATCH --output=Logs/14-macs-peaks/slurm-%A_%a.out
#SBATCH --error=Logs/14-macs-peaks/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
macs3 --version
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
# "mm" is the compiled-in Mus Musculus genome size.  See MACS3 docs.
# BAMPE = BAM Paired End
# Would it make any difference to use the merged BAMs instead?
set -x
macs3 callpeak --nomodel \
    -t Results/10-remove-duplicates/$cell_type-*-mapq1.bam \
    -f BAMPE -g mm --call-summits -n ATAC-$cell_type \
    --keep-dup all --outdir Results/14-macs-peaks
