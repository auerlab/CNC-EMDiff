#!/usr/bin/env bash

## Make peaks from PE ATAC-seq aligned (BWA) reads using MACS2
## Then we can quantify peaks using featureCounts

# Check this, just copied from Andrea's
#SBATCH --mem-per-cpu=6g
#SBATCH --output=7-macs-peaklets/slurm-%A_%a.out
#SBATCH --error=7-macs-peaklets/slurm-%A_%a.err

echo Peak calling, merging replicates across all time points, calling peaklets...
# Does order of the filenames matter?
macs2 callpeak --nomodel -t 4-bwa-mem/*_L001-nodup-uniq.bam \
    -f BAMPE -g 1.37e+09 --call-summits -n ATAC.nodup.unique.macs.peaklets \
    --keep-dup all --outdir 7-macs-peaklets
