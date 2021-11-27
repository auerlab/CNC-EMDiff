#!/bin/sh -e

#############################################################################
#   Merge BAM files for the same cell line and time point, e.g.
#
#   This means every third sample, e.g.
#       chondro-sample1-*-time1-*
#       chondro-sample4-*-time1-*
#       chondro-sample7-*-time1-*
#
#   ls Data/01-trim
#
#   One job to merge each group of 3 input bams
#   5-kallisto-quant places sample N in directory Data/05-kallisto-quant/N,
#   so one job for samples 1,4,7, another for 2,5,8, etc.
#############################################################################

#SBATCH --array=1-3,10-12
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
# top showed a steady 25M virtual, 12 res with a spike to 42/30 at the
# end (samtools index?) on FreeBSD 13
#SBATCH --mem=50m
#SBATCH --output=Logs/06-merge-kallisto-bams/slurm-%A_%a.out
#SBATCH --error=Logs/06-merge-kallisto-bams/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}

# Document software versions used for publication
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
    uname -a > Logs/06-merge-kallisto-bams/os-version-$SLURM_JOB_ID.txt 2>&1
    samtools --version > Logs/06-merge-kallisto-bams/samtools-version-$SLURM_JOB_ID.txt 2>&1
fi
hostname

dir=Data/05-kallisto-quant
input_bam1=$(echo $dir/*-sample${SLURM_ARRAY_TASK_ID}-*/pseudoalignments.bam)
input_bam2=$(echo $dir/*-sample$((SLURM_ARRAY_TASK_ID + 3))-*/pseudoalignments.bam)
input_bam3=$(echo $dir/*-sample$((SLURM_ARRAY_TASK_ID + 6))-*/pseudoalignments.bam)

# Task ID is either 1-3 (chondro) or 10-12 (neuro)
first=$SLURM_ARRAY_TASK_ID
last=$(($SLURM_ARRAY_TASK_ID + 2))
if [ $SLURM_ARRAY_TASK_ID -le 3 ]; then
    output_bam=Data/06-merge-kallisto-bams/chondro-$first-$last.bam
else
    output_bam=Data/06-merge-kallisto-bams/neuro-$first-$last.bam
fi

set -x
samtools merge -o $output_bam $input_bam1 $input_bam2 $input_bam3
samtools index $output_bam
