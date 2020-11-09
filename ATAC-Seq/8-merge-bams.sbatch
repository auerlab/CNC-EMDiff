#!/bin/sh -e

##########################################################################
#   Script description:
#       Merge alignments of replicate samples for viewing?
#
#   Usage:
#       SLURM Cluster:
#           sbatch 8-merge-bams.sbatch
#       No cluster:
#           ../slurm-sim 8-merge-bams.sbatch |& tee 1.log
#
#   History:
#       Based on the work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#   Date        Name        Modification
#   2020-03-15  Jason Bacon Begin
##########################################################################

##########################################################################
#   Main
##########################################################################

# SLURM parameters
# About 10M on Albacore
#SBATCH --mem=20
# Each job uses sample X, X+3, and X+6
#SBATCH --array=1-6
#SBATCH --output=8-merge-bams/slurm-%A_%a.out
#SBATCH --error=8-merge-bams/slurm-%A_%a.err

# : ${SLURM_ARRAY_TASK_ID:=6}

if [ -z $SLURM_ARRAY_TASK_ID ]; then
    printf "$0 must be run under sbatch or ../slurm-sim\n"
    exit 1
fi

# Document software versions used for publication
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
    uname -a > 8-merge-bams/os-version.txt 2>&1
    samtools --version > 8-merge-bams/samtools-version.txt 2>&1
fi

# Samples 1-9 are Chondrocyte and 10-18 are Neural
# How things are supposed to be:
#   Every third sample is from the same time point, e.g.
#   The # following _S is the sample ID
#   A, B, and C are time points, see:
#   Illumina Sequencing Samples_UDVADIA_OUR REFERENCE KEY, updated 082919.xlsx
#   1, 4, and 7 are Chondrocyte time 1
#   2, 5, and 8 are Chondrocyte time 2
#   10, 13, and 16 are Neural time 1
#
#   This is wrong due to a mixup in prep, time points and replicates
#   are reversed.  So for ATAC (not for RNA) A, B, and C are replicates
#   and 1, 2, and 3 are time points.
#
# CCA1A_S1_L001-nodup-uniq.bam      Replicate A, time 1
# CCA1B_S2_L001-nodup-uniq.bam      Replicate B, time 1
# CCA1C_S3_L001-nodup-uniq.bam      ...
# CCA2A_S4_L001-nodup-uniq.bam
# CCA2B_S5_L001-nodup-uniq.bam
# CCA2C_S6_L001-nodup-uniq.bam
# CCA3A_S7_L001-nodup-uniq.bam
# CCA3B_S8_L001-nodup-uniq.bam
# CCA3C_S9_L001-nodup-uniq.bam
# NCA1A_S10_L001-nodup-uniq.bam
# NCA1B_S11_L001-nodup-uniq.bam
# NCA1C_S12_L001-nodup-uniq.bam
# NCA2A_S13_L001-nodup-uniq.bam
# NCA2B_S14_L001-nodup-uniq.bam
# NCA2C_S15_L001-nodup-uniq.bam
# NCA3A_S16_L001-nodup-uniq.bam
# NCA3B_S17_L001-nodup-uniq.bam
# NCA3C_S18_L001-nodup-uniq.bam
#
# Merge BAM files for the same cell line and time point
# Example input file: CCA1A_S1_L001-nodup-uniq.bam
# Wildcards are not expanded in the RHS of a normal variable assignment
# so run a command in $()

case $SLURM_ARRAY_TASK_ID in
1)
    # Should be this, but time point and replicate were switched on samples
    # input_bams=$(ls 6-remove-duplicates/CCA*A*.bam)
    sample_group=CCA1
    ;;
2)
    # Should be this, but time point and replicate were switched on samples
    # input_bams=$(ls 6-remove-duplicates/CCA*B*.bam)
    sample_group=CCA2
    ;;
3)
    # Should be this, but time point and replicate were switched on samples
    # input_bams=$(ls 6-remove-duplicates/CCA*C*.bam)
    sample_group=CCA3
    ;;
4)
    # Should be this, but time point and replicate were switched on samples
    # input_bams=$(ls 6-remove-duplicates/NCA*A*.bam)
    sample_group=NCA1
    ;;
5)
    # Should be this, but time point and replicate were switched on samples
    # input_bams=$(ls 6-remove-duplicates/NCA*B*.bam)
    sample_group=NCA2
    ;;
6)
    # Should be this, but time point and replicate were switched on samples
    # input_bams=$(ls 6-remove-duplicates/NCA*C*.bam)
    sample_group=NCA3
    ;;
esac
time=${sample_group#*CA}
input_bams=$(ls 6-remove-duplicates/${sample_group}*.bam)

if [ $SLURM_ARRAY_TASK_ID -le 3 ]; then
    output_bam=8-merge-bams/atac-chondrocyte-$time.bam
else
    output_bam=8-merge-bams/atac-neural-$time.bam
fi

set -x
samtools merge $output_bam $input_bams
samtools index $output_bam
