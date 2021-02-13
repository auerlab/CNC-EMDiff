#!/bin/sh -e

# Go easy on the I/O.  Jobs don't take long if they're not I/O blocked.
#SBATCH --array=1-18%6 --ntasks=2
# fgrep with thousands of search strings needs a fair amount of RAM
#SBATCH --mem=1g
#SBATCH --output=6-remove-duplicates/uniquely-mapped-%A_%a.out
#SBATCH --error=6-remove-duplicates/uniquely-mapped-%A_%a.err

: ${SLURM_ARRAY_TASK_ID:=1}

cd 6-remove-duplicates
bam_files=$(ls *.bam)
my_sam_file=$(echo $bam_files | awk -v id=$SLURM_ARRAY_TASK_ID '{ print $id }')

printf "===\n$my_sam_file:\n"
printf "Total alignments: "
total=$(samtools view -c $my_sam_file)
printf "$total\n"

printf "Filtering for mapped and properly paired reads: "
# man samtools-flagstat
# properly paired: both 0x1 and 0x2 bits set and 0x4 bit not set
samtools view -b -f 0x3 -F 0x4 -o $my_sam_file.properly-paired $my_sam_file
properly_paired=$(samtools view -c $my_sam_file.properly-paired)
printf "Properly paired: $properly_paired\n"

printf "Gathering secondary read QNAMEs: "
samtools view -f 0x100 $my_sam_file.properly-paired \
    | mawk '{ print $1 }' | sort | uniq > $my_sam_file.secondary-qnames
cat $my_sam_file.secondary-qnames | wc -l

printf "Alignments excluding unmapped and multiply mapped: "
# Maybe sort by qname and step through to avoid grep memory consumption
uniquely_mapped=$(samtools view $my_sam_file.properly-paired \
    | mawk '{ print $1 }' | fgrep -v -f $my_sam_file.secondary-qnames | wc -l)
printf "$uniquely_mapped\n"

percent=$(printf "%u / %u * 100\nquit\n" $uniquely_mapped $properly_paired | bc -l)
printf "%s\t%s\t%s\t%s\t%s\n" File Total Paired UM '%'
printf "%s\t%u\t%u\t%u\t%f\n" $my_sam_file $total $properly_paired \
    $uniquely_mapped $percent | tee $my_sam_file.uniquely-mapped
