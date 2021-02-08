#!/bin/sh -e

#SBATCH --array=1-18
#SBATCH --output=6-remove-duplicates/uniquely-mapped-%A_%a.out
#SBATCH --error=6-remove-duplicates/uniquely-mapped-%A_%a.err

: ${SLURM_ARRAY_TASK_ID:=1}

cd 6-remove-duplicates
bam_files=$(ls *.bam)
my_bam_file=$(echo $bam_files | awk -v id=$SLURM_ARRAY_TASK_ID '{ print $id }')
mapped=$(samtools view -F 0x4 $my_bam_file | wc -l)
uniquely_mapped=$(samtools view -F 0x104 $my_bam_file | wc -l)
percent=$(printf "%u / %u * 100\nquit\n" $uniquely_mapped $mapped | bc -l)
printf "%s\t%u\t%u\t%f\n" $my_bam_file $mapped $uniquely_mapped $percent
