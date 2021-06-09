#!/bin/sh -e

# Samples 1-9 are Chondrocyte and 10-18 are Neural
# Every third sample is from the same time point, e.g.
# 1, 4, and 7 are Chondrocyte time 1
# 2, 5, and 8 are Chondrocyte time 2
# ...
#
# CE1A_S1_L002-R1.fastq.gz
# CE1A_S1_L002-R2.fastq.gz
# CE1B_S2_L002-R1.fastq.gz
# CE1B_S2_L002-R2.fastq.gz
# CE1C_S3_L002-R1.fastq.gz
# CE1C_S3_L002-R2.fastq.gz
# CE2A_S4_L002-R1.fastq.gz
# CE2A_S4_L002-R2.fastq.gz
# CE2B_S5_L002-R1.fastq.gz
# CE2B_S5_L002-R2.fastq.gz
# CE2C_S6_L002-R1.fastq.gz
# CE2C_S6_L002-R2.fastq.gz
#
# Merge BAM files for the same cell line and time point

# One job for input bam X, X+3, and X+6
#SBATCH --array=1-3,10-12
# FIXME: Verify memory use
#SBATCH --mem=1g
#SBATCH --output=Data/6-merge-bams/slurm-%A_%a.out
#SBATCH --error=Data/6-merge-bams/slurm-%A_%a.err

input_bam1=4-kallisto-quant/$SLURM_ARRAY_TASK_ID/pseudoalignments.bam
input_bam2=4-kallisto-quant/$((SLURM_ARRAY_TASK_ID + 3))/pseudoalignments.bam
input_bam3=4-kallisto-quant/$((SLURM_ARRAY_TASK_ID + 6))/pseudoalignments.bam

if [ $SLURM_ARRAY_TASK_ID -le 3 ]; then
    output_bam=Data/6-merge-bams/rna-chondrocyte-$SLURM_ARRAY_TASK_ID.bam
else
    output_bam=Data/6-merge-bams/rna-neural-$((SLURM_ARRAY_TASK_ID - 9)).bam
fi

samtools merge $output_bam $input_bam1 $input_bam2 $input_bam3
samtools index $output_bam
