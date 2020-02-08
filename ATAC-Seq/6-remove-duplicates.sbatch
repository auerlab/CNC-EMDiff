#!/bin/sh -e

# This script is I/O-intensive, so use slow nodes and limit concurrency
#SBATCH --array=1-18%9
# Totally lucky guess, confirmed on Peregrine
#SBATCH --mem=1g
#SBATCH --output=6-remove-duplicates/slurm-%A_%a.out
#SBATCH --error=6-remove-duplicates/slurm-%A_%a.err

# Spoof SLURM_ARRAY_TASK_ID if not set for testing outside SLURM env
: ${SLURM_ARRAY_TASK_ID:=1}

cd 6-remove-duplicates

# Variable assignments don't expand globbing patterns, so insert an ls command
saminput=$(ls ../4-bwa-mem/*_S${SLURM_ARRAY_TASK_ID}_L001.sam)
sorted=${saminput%.sam}-sorted.bam
nodup=${saminput%.sam}-nodup.bam
nodupunique=${saminput%.sam}-nodup-uniq.bam

printf 'Convert SAM to BAM, sort and index...\n'

# Clean up from prior interrupted sort processes
rm -f $sorted-*tmp*

# Sort by leftmost coordinate
samtools sort $saminput -o $sorted

samtools index $sorted
# Are these used?
samtools idxstats $sorted > $sorted.idxstats
samtools flagstat $sorted > $sorted.flagstat

# Can we remove duplicates before aligning?
# Seems like this might save a lot of time on redundant alignments
printf 'Removing duplicate reads...\n'

# fixmate requires name-sorted input
samtools sort -n -o $sorted-namesort.bam $sorted 
rm -f $sorted $sorted.*

samtools fixmate -m $sorted-namesort.bam $sorted-fixmate.bam
rm -f $sorted-namesort.bam

# markdup requires coordinate-sorted input
samtools sort -o $sorted-fixmate-sort.bam $sorted-fixmate.bam
samtools markdup -l 100 -r -s $sorted-fixmate-sort.bam $nodup
rm -f $sorted-fixmate-sort.bam

samtools index $nodup
# Are these used?
samtools idxstats $nodup > $nodup.idxstats
samtools flagstat $nodup > $nodup.flagstat

printf 'Removing multiple mapped reads...\n'
# Default quality min (-q) is 0
samtools view -bq 1 $nodup > $nodupunique
rm -f $nodup $nodup.*

samtools index $nodupunique
samtools idxstats $nodupunique > $nodupunique.idxstats
samtools flagstat $nodupunique > $nodupunique.flagstat
