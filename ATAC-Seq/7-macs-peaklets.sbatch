#!/usr/bin/env bash

##########################################################################
#   Script description:
#       Call peaks in alignments PE ATAC-seq aligned (BWA) reads
#
#   Usage:
#       SLURM Cluster:
#           sbatch 7-macs-peaklets.sbatch
#       No cluster:
#           ../slurm-sim 7-macs-peaklets.sbatch 1 2 |& tee 1.log
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
#SBATCH --output=7-macs-peaklets/slurm-%A_%a.out
#SBATCH --error=7-macs-peaklets/slurm-%A_%a.err

if [ -z $SLURM_ARRAY_TASK_ID ]; then
    printf "$0 must be run under sbatch or ../slurm-sim\n"
    exit 1
fi

# Document software versions used for publication
uname -a > 7-macs-peaks/os-version.txt 2>&1
macs2 --version > 7-macs-peaks/macs2-version.txt 2>&1

case $task_id in
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

printf "$cell_type: Peak calling, merging replicates across all time points, calling peaklets...\n"
# "mm" is the compiled-in Mus Musculus genome size.  See MACS2 docs.
# BAMPE = BAM Paired End
# Does order of the filenames matter?
# Could we use the merged BAMs from 8- instead?
macs2 callpeak --nomodel -t 4-bwa-mem/$cell_type*_L001-nodup-uniq.bam \
    -f BAMPE -g mm --call-summits -n ATAC-$cell_type-macs.peaklets \
    --keep-dup all --outdir 7-macs-peaklets
