#!/bin/sh -e

##########################################################################
#   Script description:
#       Trim adapters and low quality ends from Lumina reads
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#   History:
#   Date        Name        Modification
#   2019-09-11  Jason Bacon Begin
##########################################################################

##########################################################################
#   Main
##########################################################################

# Each job in the array should get 2 cores on the same node to accomodate
# the cutadapt process and companion pigz/gzip process
#SBATCH --array=1-18 --cpus-per-task=2
#SBATCH --mem-per-cpu=100
#SBATCH --output=Data/1-trim/slurm-%A_%a.out
#SBATCH --error=Data/1-trim/slurm-%A_%a.err

: ${SLURM_ARRAY_TASK_ID:=1}

input1=$(ls ../Raw/190822_AHFN3KDRXX-Lane2/*_S${SLURM_ARRAY_TASK_ID}_*R1*)
input2=$(ls ../Raw/190822_AHFN3KDRXX-Lane2/*_S${SLURM_ARRAY_TASK_ID}_*R2*)
base_name=$(basename $input1)
stem=${base_name%%_R*.fastq.xz}

# TrimGalore is just a wrapper around cutadapt and it cannot cut a fixed
# number of bases at the same time as trimming adapters ( -u and -U in
# cutadapt ), so we run cutadapt directly.

#trim_galore --illumina --stringency 3 -q 20 --paired -o Data/1-trim \
#    $input1 $input2 > Data/1-trim/${stem}-trash.out
#trim_galore --illumina --stringency 3 -q 20 --paired --trim1 -o Data/1-trim-trim1 \
#    $input1 $input2 > Data/1-trim-trim1/${stem}-trash.out

# Command from trim_galore log:
#   cutadapt -j 1 -e 0.1 -q 20 -O 3 -a AGATCGGAAGAGC
#   -j 1 -O 3 -e 0.1 are default for cutadapt

# Remove bases at 5' end due to bias and last base at 3' end
# 5' bias is probably due to non-randomness of start sequences and not a
# quality issue, so -u +15 -U +15 is probably not necessary
# -u +15 -U +15 \

cutadapt \
    --quality-cutoff=24 \
    -u -3 -U -3 \
    -a AGATCGGAAGAGCACAC -A AGATCGGAAGAGCGTCG \
    --minimum-length=30 \
    --output=Data/1-trim/${stem}-R1.fastq.xz \
    --paired-output=Data/1-trim/${stem}-R2.fastq.xz \
    $input1 $input2
