#!/bin/sh -e

##########################################################################
#   Script description:
#       Filter out duplicate reads
#       This script is I/O-intensive, so use slow nodes and limit
#       concurrency.
#
#   Usage:
#       SLURM Cluster:
#           sbatch 6-remove-duplicates.sbatch
#       No cluster:
#           ./slurm-sim 6-remove-duplicates.sbatch 1 18 |& tee 6.log
#
#   History:
#       Based on the work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#   Date        Name        Modification
#   2020-03-01  Jason Bacon Begin
##########################################################################

##########################################################################
#   Main
##########################################################################

# SLURM parameters
#SBATCH --array=1-18%9
#SBATCH --mem=2g
#SBATCH --output=6-remove-duplicates/slurm-%A_%a.out
#SBATCH --error=6-remove-duplicates/slurm-%A_%a.out

if [ -z $SLURM_ARRAY_TASK_ID ]; then
    printf "$0 must be run under sbatch or ./slurm-sim\n"
    exit 1
fi

# Document software versions used for publication
uname -a > 6-remove-duplicates/os-version.txt 2>&1
samtools --version > 6-remove-duplicates/samtools-version.txt 2>&1

# FIXME: Remove duplicate reads before alignment to cut down on redundant
# alignments.  Still need to run this step since we cannot identify all
# duplicate alignments before alignment.

cd 6-remove-duplicates

# Variable assignments don't expand globbing patterns, so insert an ls command
saminput=$(ls ../4-bwa-mem/*_S${SLURM_ARRAY_TASK_ID}_L001.sam)
sorted=${saminput%.sam}-sorted.bam
nodup=${saminput%.sam}-nodup.bam
nodupunique=${saminput%.sam}-nodup-uniq.bam
printf "$saminput $sorted $nodup $nodupunique\n"

# Clean up from prior interrupted sort processes
rm -f $sorted-*tmp*

# Is there a reason to sort by coordinates before sorting by name?
# Sort by leftmost coordinate
# samtools sort $saminput -o $sorted
# samtools index $sorted
# Are these even used?
# samtools idxstats $sorted > $sorted.idxstats
# samtools flagstat $sorted > $sorted.flagstat

# fixmate requires name-sorted input
# Apparently must be in a file, so we can't pipe sorted input to it?
# samtools sort -n -o $sorted-namesort.bam $sorted 
# -n: Sort by QNAME
# -m: Specify max memory use.  More mem means fewer temp files
printf 'Sorting by name for fixmate...\n'
samtools sort -n -m 2g -o $sorted-namesort.bam $saminput
rm -f $sorted $sorted.*

# -m: Add mate score tags to help markdup select the best reads to keep
printf 'Adding mate scores for markdup with fixmate...\n'
samtools fixmate -m $sorted-namesort.bam $sorted-fixmate.bam
rm -f $sorted-namesort.bam

# markdup requires coordinate-sorted input
# -m: Specify max memory use.  More mem means fewer temp files
printf 'Sorting by coordinates for markdup...\n'
samtools sort -m 2g -o $sorted-fixmate-sort.bam $sorted-fixmate.bam
rm -f $sorted-fixmate.bam

# -l: Expected read length
# -r: Remove duplicates
# -s: Print basic stats
printf 'Removing duplicate reads with markdup...\n'
samtools markdup -l 100 -r -s $sorted-fixmate-sort.bam $nodup
rm -f $sorted-fixmate-sort.bam

printf 'Generating index and stats...\n'
# Are these used?
samtools index $nodup
samtools idxstats $nodup > $nodup.idxstats
samtools flagstat $nodup > $nodup.flagstat

printf 'Removing reads with MAPQ < 1...\n'
# Default quality min (-q) is 0
# -b: Output BAM format
# -q: Minimum MAPQ value
samtools view -b -q 1 $nodup > $nodupunique
rm -f $nodup $nodup.*

printf 'Generating index and stats...\n'
samtools index $nodupunique
samtools idxstats $nodupunique > $nodupunique.idxstats
samtools flagstat $nodupunique > $nodupunique.flagstat
