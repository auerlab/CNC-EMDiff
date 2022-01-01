#!/bin/sh -e

#############################################################################
#   Description:
#       Merge BAM files for the same cell line and time point, e.g.
#
#       This means every third sample, e.g.
#           chondro-sample1-*-time1-*
#           chondro-sample4-*-time1-*
#           chondro-sample7-*-time1-*
#
#       ls Data/03-trim
#
#       One job to merge each group of 3 input bams
#       so one job for samples 1,4,7, another for 2,5,8, etc.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#############################################################################

# Each task will merge 3 BAMs from samples n, n+3, and n+6
#SBATCH --array=1-3,10-12
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
# top showed a steady 25M virtual, 12 res with a spike to 42/30 at the
# end (samtools index?) on FreeBSD 13
#SBATCH --mem=50m
#SBATCH --output=Logs/17-merge-hisat2-bams/slurm-%A_%a.out
#SBATCH --error=Logs/17-merge-hisat2-bams/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
samtools --version
pwd

dir=Data/16-hisat2-align
input_bam1=$(echo $dir/*-sample${SLURM_ARRAY_TASK_ID}-*.bam)
input_bam2=$(echo $dir/*-sample$((SLURM_ARRAY_TASK_ID + 3))-*.bam)
input_bam3=$(echo $dir/*-sample$((SLURM_ARRAY_TASK_ID + 6))-*.bam)

# Task ID is either 1-3 (chondro) or 10-12 (neuro)
first=$SLURM_ARRAY_TASK_ID
last=$(($SLURM_ARRAY_TASK_ID + 2))
if [ $SLURM_ARRAY_TASK_ID -le 3 ]; then
    output_bam=Data/17-merge-hisat2-bams/chondro-$first-$last.bam
else
    output_bam=Data/17-merge-hisat2-bams/neuro-$first-$last.bam
fi

set -x
samtools merge -o $output_bam $input_bam1 $input_bam2 $input_bam3
samtools index $output_bam
