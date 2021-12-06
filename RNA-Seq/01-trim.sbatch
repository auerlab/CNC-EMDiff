#!/bin/sh -e

##########################################################################
#   Script description:
#       Trim adapters and low quality ends from Lumina reads
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#
#   History:
#   Date        Name        Modification
#   2019-09-11  Jason Bacon Begin
##########################################################################

# Each job in the array will run a cutadapt (python) process and a
# compression process for part of the time.  If you don't want to
# oversubscribe compute nodes even for a little while, add --cpus-per-task=2
#SBATCH --array=1-18
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem-per-cpu=250
#SBATCH --output=Logs/01-trim/slurm-%A_%a.out
#SBATCH --error=Logs/01-trim/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
cutadapt --version
pwd

input1=$(ls Data/Raw-renamed/*-sample${SLURM_ARRAY_TASK_ID}-*R1*)
input2=$(ls Data/Raw-renamed/*-sample${SLURM_ARRAY_TASK_ID}-*R2*)
base=$(basename $input1)
stem=${base%%-R*.fastq.xz}

##############################################################################
# TrimGalore is just a wrapper around cutadapt and it cannot cut a fixed
# number of bases at the same time as trimming adapters ( -u and -U in
# cutadapt ), so we run cutadapt directly.
#
# trim_galore --illumina --stringency 3 -q 20 --paired -o Data/01-trim \
#    $input1 $input2 > Data/01-trim/${stem}-trash.out
# trim_galore --illumina --stringency 3 -q 20 --paired --trim1 -o Data/01-trim-trim1 \
#    $input1 $input2 > Data/01-trim-trim1/${stem}-trash.out
#
# Command from trim_galore log:
#   cutadapt -j 1 -e 0.1 -q 20 -O 3 -a AGATCGGAAGAGC
#   -j 1 -O 3 -e 0.1 are default for cutadapt

##############################################################################
# Remove ~15 bases at 5' end due to bias and last base at 3' end
# 5' bias is probably due to non-random cleavage and not a
# quality issue, so -u +15 -U +15 is probably not necessary
# -u +15 -U +15 \

cutadapt \
    --quality-cutoff=24 \
    -u -3 -U -3 \
    -a AGATCGGAAGAGCACAC -A AGATCGGAAGAGCGTCG \
    --minimum-length=30 \
    --output=Data/01-trim/${stem}-R1.fastq.xz \
    --paired-output=Data/01-trim/${stem}-R2.fastq.xz \
    $input1 $input2
