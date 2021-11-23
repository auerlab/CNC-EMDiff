#!/bin/sh -e

#############################################################################
#   Merge BAM files for the same cell line and time point, e.g.
#
#   chondro-sample*-rep*-time1-R1.fastq.xz
#
#   One job for input bam X, X+3, and X+6
#   5-kallisto-quant places sample N in directory Data/5-kallisto-quant/N
#############################################################################

#SBATCH --array=1-3,10-12
# FIXME: Verify memory use
#SBATCH --mem=1g
#SBATCH --output=Data/6-merge-bams/slurm-%A_%a.out
#SBATCH --error=Data/6-merge-bams/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}

# Document software versions used for publication
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
    uname -a > Data/6-merge-bams/os-version-$SLURM_JOB_ID.txt 2>&1
    samtools --version > Data/6-merge-bams/samtools-version-$SLURM_JOB_ID.txt 2>&1
fi

input_bam1=Data/5-kallisto-quant/$SLURM_ARRAY_TASK_ID/pseudoalignments.bam
input_bam2=Data/5-kallisto-quant/$((SLURM_ARRAY_TASK_ID + 3))/pseudoalignments.bam
input_bam3=Data/5-kallisto-quant/$((SLURM_ARRAY_TASK_ID + 6))/pseudoalignments.bam

# Task ID is either 1-3 (chondro) or 10-12 (neuro)
if [ $SLURM_ARRAY_TASK_ID -le 3 ]; then
    output_bam=Data/6-merge-bams/rna-chondrocyte-$SLURM_ARRAY_TASK_ID.bam
else
    output_bam=Data/6-merge-bams/rna-neural-$((SLURM_ARRAY_TASK_ID - 9)).bam
fi

samtools merge $output_bam $input_bam1 $input_bam2 $input_bam3
samtools index $output_bam
