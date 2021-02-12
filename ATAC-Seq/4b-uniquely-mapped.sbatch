#!/bin/sh -e

#SBATCH --array=1-18
#SBATCH --output=4-bwa-mem/uniquely-mapped-%A_%a.out
#SBATCH --error=4-bwa-mem/uniquely-mapped-%A_%a.err

: ${SLURM_ARRAY_TASK_ID:=1}

cd 4-bwa-mem
sam_files=$(ls *.sam)
my_sam_file=$(echo $sam_files | awk -v id=$SLURM_ARRAY_TASK_ID '{ print $id }')

# Use this instead of total?
# mapped=$(samtools view -F 0x4 $my_sam_file | wc -l)

printf "===\n$my_sam_file:\n"
printf "Total alignments: "
total=$($sam_reader $my_sam_file | mawk '{ print $1 }' | wc -l)
printf "$total\n"

printf "Gathering unmapped read QNAMEs...\n"
samtools view -f 0x4 $my_sam_file | mawk '{ print $1 }' | sort | uniq > unmapped-qnames.txt
wc -l unmapped-qnames.txt

printf "Gathering secondary read QNAMEs...\n"
samtools view -f 0x100 $my_sam_file | mawk '{ print $1 }' | sort | uniq > secondary-qnames.txt
wc -l secondary-qnames.txt

printf "Alignments excluding unmapped and multiply mapped: "
cat unmapped-qnames.txt secondary-qnames.txt > remove.txt
uniquely_mapped=$($sam_reader $my_sam_file | mawk '{ print $1 }' | fgrep -v -f remove.txt | wc -l)

percent=$(printf "%u / %u * 100\nquit\n" $uniquely_mapped $total | bc -l)
printf "%s\t%u\t%u\t%f\n" $my_sam_file $total $uniquely_mapped $percent
