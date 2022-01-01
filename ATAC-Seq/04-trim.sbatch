#!/bin/sh -e

##########################################################################
#   Script description:
#       Trim adapters and low quality ends from Lumina reads
#
#   Usage:
#       SLURM Cluster:
#           sbatch 1-trim.sbatch
#       No cluster:
#           ../slurm-sim 1-trim.sbatch |& tee 1.log
#
#   History:
#       Based on the work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#   Date        Name        Modification
#   2020-02-02  Jason Bacon Begin
##########################################################################

##########################################################################
#   Main
##########################################################################

# 2 cores per task to accommodate pigz
# We don't get full utilization of both cores, but this speeds up the
# job somewhat.  If you prefer to maximize CPU utilization, set
# cpus-per-task to 1.
#SBATCH --array=1-18 --cpus-per-task=2
#SBATCH --mem-per-cpu=300
#SBATCH --output=Logs/01-trim/slurm-%A_%a.out
#SBATCH --error=Logs/01-trim/slurm-%A_%a.err

# Defaults if not set by SLURM or user
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_CPUS_PER_TASK:=2}

# Document software versions used for publication
uname -a
cutadapt --version

# One iteration if running under SLURM, all iterations otherwise
printf "Running task $SLURM_ARRAY_TASK_ID...\n"

# chondro-sample1-rep1-time1-R1.fastq.xz
input1="Data/Raw-renamed/*-sample${SLURM_ARRAY_TASK_ID}-*R1*"
input2="Data/Raw-renamed/*-sample${SLURM_ARRAY_TASK_ID}-*R2*"
base_name=$(basename $input1)
stem=${base_name%%-R*.fastq.xz}

# TrimGalore cannot cut a fixed number of bases while trimming adapters.
# We want to remove bases at 5' end due to bias and last base at 3' end

# Andrea's code:
# trim_galore --stringency 3 -q 20 --paired -o 1-trim \
#    $input1 $input2 > 1-trim/${stem}-trash.out
# trim_galore --stringency 3 -q 20 --paired --trim1 -o 1-trim-trim1 \
#    $input1 $input2 > 1-trim-trim1/${stem}-trash.out

# From trim_galore log:
#   cutadapt -j 1 -e 0.1 -q 20 -O 3 -a AGATCGGAAGAGC
# -j 1 -O 3 -e 0.1 are defaults for cutadapt

# ATAC-Seq adapter: https://www.biostars.org/p/367067/
adapter=CTGTCTCTTATACACATCT

# Strong bias in first ~15 bases
# Most quality scores are above 30, so we can be pickier thant the
# traditional cut-off of 20
# Not sure about min length for ATAC
# Need 2 cores for each cutadapt process to accommodate pigz/xzip

set -x
cutadapt \
    --cores=$SLURM_CPUS_PER_TASK \
    --quality-cutoff=24 \
    -u -3 -U -3 \
    -a $adapter -A $adapter \
    --minimum-length=20 \
    --output=Data/01-trim/${stem}-R1.fastq.xz \
    --paired-output=Data/01-trim/${stem}-R2.fastq.xz \
    $input1 $input2
