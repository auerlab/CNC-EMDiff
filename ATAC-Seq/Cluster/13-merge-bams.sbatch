#!/bin/sh -e

##########################################################################
#   Script description:
#       Merge replicates (BAM files for the same cell line and time point)
#       mainly for viewing in IGV
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

# About 10M on Albacore
#SBATCH --mem=20
# 18 samples, merging 3 time points
#SBATCH --array=1-6
#SBATCH --output=Logs/13-merge-bams/slurm-%A_%a.out
#SBATCH --error=Logs/13-merge-bams/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
samtools --version
pwd

case $SLURM_ARRAY_TASK_ID in
1)
    cell_line=chondro
    time=1
    ;;
2)
    cell_line=chondro
    time=2
    ;;
3)
    cell_line=chondro
    time=3
    ;;
4)
    cell_line=neuro
    time=1
    ;;
5)
    cell_line=neuro
    time=2
    ;;
6)
    cell_line=neuro
    time=3
    ;;
esac
input_bams=$(ls Results/10-remove-duplicates/$cell_line-*-time$time-*-mapq1.bam)

if [ $SLURM_ARRAY_TASK_ID -le 3 ]; then
    output_bam=Results/13-merge-bams/atac-chondro-$time.bam
else
    output_bam=Results/13-merge-bams/atac-neuro-$time.bam
fi

set -x
samtools merge -f $output_bam $input_bams
samtools index $output_bam
