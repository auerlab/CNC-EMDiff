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
#SBATCH --output=1-trim/slurm-%A_%a.out
#SBATCH --error=1-trim/slurm-%A_%a.err

if [ -z $SLURM_ARRAY_TASK_ID ]; then
    printf "$0 must be run under sbatch or ../slurm-sim\n"
    exit 1
fi

# Defaults if not set by SLURM or user
: ${SLURM_CPUS_PER_TASK:=2}

# Document software versions used for publication
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
    uname -a > 1-trim/os-version.txt 2>&1
    cutadapt --version > 1-trim/cutadapt-version.txt 2>&1
fi

# One iteration if running under SLURM, all iterations otherwise
printf "Running task $SLURM_ARRAY_TASK_ID...\n"

input1="../Raw/190822_AHFN3KDRXX-Lane1/*_S${SLURM_ARRAY_TASK_ID}_*R1*"
input2="../Raw/190822_AHFN3KDRXX-Lane1/*_S${SLURM_ARRAY_TASK_ID}_*R2*"
base_name=$(basename $input1)
stem=${base_name%%_R*.fastq.xz}

# TrimGalore cannot cut a fixed number of bases while trimming adapters.
# We want to remove bases at 5' end due to bias and last base at 3' end

# Andrea's code:
# trim_galore --illumina --stringency 3 -q 20 --paired -o 1-trim \
#    $input1 $input2 > 1-trim/${stem}-trash.out
# trim_galore --illumina --stringency 3 -q 20 --paired --trim1 -o 1-trim-trim1 \
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


cutadapt \
    --cores=$SLURM_CPUS_PER_TASK \
    --quality-cutoff=24 \
    -u -3 -U -3 \
    -u +15 -U +15 \
    -a $adapter -A $adapter \
    --minimum-length=20 \
    --output=1-trim/${stem}-R1.fastq.xz \
    --paired-output=1-trim/${stem}-R2.fastq.xz \
    $input1 $input2
