#!/bin/sh -e

# Go easy on the I/O.  Jobs don't take long if they're not I/O blocked.
#SBATCH --array=1-18%2 --ntasks=2
# fgrep with millions of search strings needs lots of RAM
#SBATCH --mem=20g
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
total=$(cat $my_sam_file | mawk '{ print $1 }' | wc -l)
printf "$total\n"

printf "Gathering unmapped read QNAMEs...\n"
samtools view -f 0x4 $my_sam_file | mawk '{ print $1 }' | sort | uniq > $my_sam_file.unmapped-qnames
wc -l $my_sam_file.unmapped-qnames

printf "Gathering secondary read QNAMEs...\n"
samtools view -f 0x100 $my_sam_file | mawk '{ print $1 }' | sort | uniq > $my_sam_file.secondary-qnames
wc -l $my_sam_file.secondary-qnames

printf "Alignments excluding unmapped and multiply mapped: "
cat $my_sam_file.unmapped-qnames $my_sam_file.secondary-qnames > $my_sam_file.remove-qnames
uniquely_mapped=$(cat $my_sam_file | mawk '{ print $1 }' | fgrep -v -f $my_sam_file.remove-qnames | wc -l)

percent=$(printf "%u / %u * 100\nquit\n" $uniquely_mapped $total | bc -l)
printf "%s\t%u\t%u\t%f\n" $my_sam_file $total $uniquely_mapped $percent \
    | tee $my_sam_file.uniquely-mapped
