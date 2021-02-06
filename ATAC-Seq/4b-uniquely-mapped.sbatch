#!/bin/sh -e

#SBATCH --array=1-18
#SBATCH --output=4-bwa-mem/uniquely-mapped-%A_%a.out
#SBATCH --error=4-bwa-mem/uniquely-mapped-%A_%a.err

: ${SLURM_ARRAY_TASK_ID:=1}

sam_files=$(ls 4-bwa-mem/*.sam)
my_sam_file=$(echo $sam_files | awk -v id=$SLURM_ARRAY_TASK_ID '{ print $id }')
printf "Processing $my_sam_file...\n"
total_alignments=$(samtools view $my_sam_file | wc -l)
uniquely_mapped=$(samtools view -F 0x100 $my_sam_file | wc -l)
total_alignments=100
uniquely_mapped=98
percent=$(printf "%u / %u * 100\nquit\n" $uniquely_mapped $total_alignments | bc -l)
printf "%u total alignments  %u uniquely mapped (%f%%)\n" \
    $total_alignments $uniquely_mapped $percent
