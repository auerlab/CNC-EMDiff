#!/bin/sh -e

# The digit in the 3rd character of the filename (1,2,3) indicates replicate
# The letter in the 4th character (A, B, C) represents time point
#   A = time 0, B = time 1, C = time 2
# The number following _S is the sample #
#
# Samples 1-9 are Chondrocyte and 10-18 are Neural
# Every third sample (those with the same letter A, B, or C in the 4th
# character of the filename) is from the same time point, e.g.
#
# CE1A_S1, CE2A_S4, and CE3A_S7 are chondrocyte time 0
# CE1B_S2, CE2B_S5, and CE3B_S8 are Chondrocyte time 1
#
# CE1A_S1_L002-R1.fastq.xz
# CE1A_S1_L002-R2.fastq.xz
# CE1B_S2_L002-R1.fastq.xz
# CE1B_S2_L002-R2.fastq.xz
# CE1C_S3_L002-R1.fastq.xz
# CE1C_S3_L002-R2.fastq.xz
# CE2A_S4_L002-R1.fastq.xz
# CE2A_S4_L002-R2.fastq.xz
# CE2B_S5_L002-R1.fastq.xz
# CE2B_S5_L002-R2.fastq.xz
# CE2C_S6_L002-R1.fastq.xz
# CE2C_S6_L002-R2.fastq.xz
# CE3A_S7_L002-R1.fastq.xz
# CE3A_S7_L002-R2.fastq.xz
# CE3B_S8_L002-R1.fastq.xz
# CE3B_S8_L002-R2.fastq.xz
# CE3C_S9_L002-R1.fastq.xz
# CE3C_S9_L002-R2.fastq.xz
# NE1A_S10_L002-R1.fastq.xz
# NE1A_S10_L002-R2.fastq.xz
# ...
# NE3C_S18_L002-R2.fastq.xz

# Merge BAM files for the same cell line and time point, e.g. CE*A*
# One job for input bam X, X+3, and X+6
# 5-kallisto-quant places sample N in directory Data/5-kallisto-quant/N

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
