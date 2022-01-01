#!/bin/sh -e

##########################################################################
#   Script description:
#       Filter out duplicate reads
#       This script is I/O-intensive, so use slow nodes and limit
#       concurrency.
#
#       FIXME: Explore using Picard MarkDuplicates for this instead
#
#   Usage:
#       SLURM Cluster:
#           sbatch 6-remove-duplicates.sbatch
#       No cluster:
#           ../slurm-sim 6-remove-duplicates.sbatch |& tee 6.log
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
#SBATCH --array=1-18%6
#SBATCH --mem=2g
#SBATCH --output=6-remove-duplicates/slurm-%A_%a.out
#SBATCH --error=6-remove-duplicates/slurm-%A_%a.out

if [ -z $SLURM_ARRAY_TASK_ID ]; then
    printf "$0 must be run under sbatch or ../slurm-sim\n"
    exit 1
fi

# Document software versions used for publication
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
    uname -a > 6-remove-duplicates/os-version.txt 2>&1
    samtools --version > 6-remove-duplicates/samtools-version.txt 2>&1
fi

# FIXME: Remove duplicate reads before alignment to cut down on redundant
# alignments.  Still need to run this step since we cannot identify all
# duplicate alignments before alignment.

cd 6-remove-duplicates

# This should match only 1 file.  '*' used only to match irrelevant variations.
# The integer following 'S' in the filename should be unique (1 - 18)
# Variable assignments don't expand globbing patterns, so insert an ls command
saminput=$(ls ../4-bwa-mem/*_S${SLURM_ARRAY_TASK_ID}_L001.sam)
sambase=$(basename $saminput)
sorted=${sambase%.sam}-sorted.bam
nodup=${sambase%.sam}-nodup.bam
nodupmapq1=${sambase%.sam}-nodup-mapq1.bam
printf "saminput=$saminput, sorted=$sorted, nodup=$nodup, nodupmapq1=$nodupmapq1\n"

# Clean up from prior interrupted sort processes
rm -f $sorted-*tmp*

# Is there a reason to sort by coordinates before sorting by name?
# Sort by leftmost coordinate
# samtools sort $saminput -o $sorted
# samtools index $sorted
# Are these even used?
# samtools idxstats $sorted > $sorted.idxstats
# samtools flagstat $sorted > $sorted.flagstat

set -x

# fixmate requires name-sorted input
# Apparently must be in a file, so we can't pipe sorted input to it?
# samtools sort -n -o $sorted-namesort.bam $sorted 
# -n: Sort by QNAME
# -m: Specify max memory use.  More mem means fewer temp files
samtools sort -n -m 2g -o $sorted-namesort.bam $saminput
rm -f $sorted $sorted.*

# -m: Add mate score tags to help markdup select the best reads to keep
samtools fixmate -m $sorted-namesort.bam $sorted-fixmate.bam
# Keep name-sorted file for further filtering
# rm -f $sorted-namesort.bam

# markdup requires coordinate-sorted input
# -m: Specify max memory use.  More mem means fewer temp files
samtools sort -m 2g -o $sorted-fixmate-sort.bam $sorted-fixmate.bam
rm -f $sorted-fixmate.bam

# Remove duplicate reads with markdup
# -l: Expected read length
# -r: Remove duplicates
# -s: Print basic stats
samtools markdup -l 100 -r -s $sorted-fixmate-sort.bam $nodup
rm -f $sorted-fixmate-sort.bam

# Generating index and stats
# Are these used?
samtools index $nodup
samtools idxstats $nodup > $nodup.idxstats
samtools flagstat $nodup > $nodup.flagstat

# Remove reads with MAPQ < 1
# Default quality min (-q) is 0
# -b: Output BAM format
# -q: Minimum MAPQ value
samtools view -b -q 1 $nodup > $nodupmapq1
rm -f $nodup $nodup.*

# Generating index and stats
samtools index $nodupmapq1
samtools idxstats $nodupmapq1 > $nodupmapq1.idxstats
samtools flagstat $nodupmapq1 > $nodupmapq1.flagstat
