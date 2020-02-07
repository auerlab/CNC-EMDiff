#!/bin/sh -e

#SBATCH --array=1-18%9
# Totally lucky guess, confirmed on Peregrine
#SBATCH --mem=1g
#SBATCH --output=6-remove-duplicates/slurm-%A_%a.out
#SBATCH --error=6-remove-duplicates/slurm-%A_%a.err

# Spoof SLURM_ARRAY_TASK_ID if not set for testing outside SLURM env
: ${SLURM_ARRAY_TASK_ID:=1}

cd 6-remove-duplicates

# FIXME: Remove .save after benchmarking BWA
# Variable assignments don't expand globbing patterns, so insert an ls command
saminput=$(ls ../4-bwa-mem.save/*_S${SLURM_ARRAY_TASK_ID}_L001.sam)
sorted=${saminput%.sam}-sorted.bam
nodup=${saminput%.sam}-nodup.bam
nodupunique=${saminput%.sam}-nodup-uniq.bam

printf 'Convert SAM to BAM, sort and index...\n'
# Sort by leftmost coordinate
samtools sort $saminput -o $sorted

samtools index $sorted
samtools idxstats $sorted > $sorted.idxstats
samtools flagstat $sorted > $sorted.flagstat

# Can we remove duplicates before aligning?
# Seems like this might save a lot of time on redundant alignments
printf 'Removing duplicate reads...\n'

# fixmate requires name-sorted input
samtools sort -n -o $sorted-namesort.bam $sorted 
samtools fixmate -m $sorted-namesort.bam $sorted-fixmate.bam

# markdup requires coordinate-sorted input
samtools sort -o $sorted-fixmate-sort.bam $sorted-fixmate.bam
samtools markdup -l 100 -r -s $sorted-fixmate-sort.bam $nodup

samtools index $nodup
samtools idxstats $nodup > $nodup.idxstats
samtools flagstat $nodup > $nodup.flagstat

printf 'Removing multiple mapped reads...\n'
# Default quality min (-q) is 0
samtools view -bq 1 $nodup > $nodupunique

samtools index $nodupunique
samtools idxstats $nodupunique > $nodupunique.idxstats
samtools flagstat $nodupunique > $nodupunique.flagstat
