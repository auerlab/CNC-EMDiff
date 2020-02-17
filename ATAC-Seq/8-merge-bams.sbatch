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

# Verify memory use
#SBATCH --mem=1g
# Each job uses sample X, X+3, and X+6
#SBATCH --array=1-3,10-12
#SBATCH --output=8-merge-bams/slurm-%A_%a.out
#SBATCH --error=8-merge-bams/slurm-%A_%a.err

# Not CPU-bound, so let it use the slower nodes
##SBATCH --exclude=compute-[001-008,012]

# Example input file: CCA1A_S1_L001-nodup-uniq.bam
# The # following _S is the sample ID
# Wildcards are not expanded in the RHS of a normal variable assignment
# so run a command in $()
input_bam1=`echo 4-bwa-mem/*_S${SLURM_ARRAY_TASK_ID}_L001-nodup-uniq.bam`
input_bam2=`echo 4-bwa-mem/*_S$((SLURM_ARRAY_TASK_ID + 3))_L001-nodup-uniq.bam`
input_bam3=`echo 4-bwa-mem/*_S$((SLURM_ARRAY_TASK_ID + 6))_L001-nodup-uniq.bam`

if [ $SLURM_ARRAY_TASK_ID -le 3 ]; then
    output_bam=8-merge-bams/atac-chondrocyte-$SLURM_ARRAY_TASK_ID.bam
else
    output_bam=8-merge-bams/atac-neural-$((SLURM_ARRAY_TASK_ID - 9)).bam
fi
samtools merge $output_bam $input_bam1 $input_bam2 $input_bam3
samtools index $output_bam
