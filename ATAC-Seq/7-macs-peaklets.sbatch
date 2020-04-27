#!/usr/bin/env bash

## Make peaks from PE ATAC-seq aligned (BWA) reads using MACS2
## Then we can quantify peaks using featureCounts

# Mem requirement observed on Peregrine
#SBATCH --array=1-2
#SBATCH --mem=4g
#SBATCH --output=7-macs-peaklets/slurm-%A_%a.out
#SBATCH --error=7-macs-peaklets/slurm-%A_%a.err

# Spoof SLURM_ARRAY_TASK_ID if not set for testing outside SLURM env
: ${SLURM_ARRAY_TASK_ID:=1}

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

echo Peak calling, merging replicates across all time points, calling peaklets...
# "mm" is the compiled-in Mus Musculus genome size.  See MACS2 docs.
# BAMPE = BAM Paired End
# Does order of the filenames matter?
# Could we use the merged BAMs from 8- instead?
macs2 callpeak --nomodel -t 4-bwa-mem/$cell_type*_L001-nodup-uniq.bam \
    -f BAMPE -g mm --call-summits -n ATAC-$cell_type-macs.peaklets \
    --keep-dup all --outdir 7-macs-peaklets
