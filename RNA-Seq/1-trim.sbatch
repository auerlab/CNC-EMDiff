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
#SBATCH --output=1-trim/slurm-%A_%a.out
#SBATCH --error=1-trim/slurm-%A_%a.err

# Mortimer pkgsrc
if [ -e /sharedapps/pkg-2019Q2/bin/cutadapt ]; then
    # sh is bash on CentOS, so this is OK
    source /etc/bashrc
    module load /sharedapps/pkg-2019Q2/etc/modulefiles/pkgsrc/2019Q2
fi

input1="../Raw/190822_AHFN3KDRXX-Lane2/*_S${SLURM_ARRAY_TASK_ID}_*R1*"
input2="../Raw/190822_AHFN3KDRXX-Lane2/*_S${SLURM_ARRAY_TASK_ID}_*R2*"
base_name=$(basename $input1)
stem=${base_name%%_R*.fastq.gz}

# TrimGalore cannot cut a fixed number of bases at the same time as trimming
# adapters.  ( -u and -U in cutadapt )
#trim_galore --illumina --stringency 3 -q 20 --paired -o 1-trim \
#    $input1 $input2 > 1-trim/${stem}-trash.out
#trim_galore --illumina --stringency 3 -q 20 --paired --trim1 -o 1-trim-trim1 \
#    $input1 $input2 > 1-trim-trim1/${stem}-trash.out

# From trim_galore log: cutadapt -j 1 -e 0.1 -q 20 -O 3 -a AGATCGGAAGAGC
# -j 1 -O 3 -e 0.1 are default for cutadapt
# Remove bases at 5' end due to bias and last base at 3' end

set -x
cutadapt \
    --quality-cutoff=24 \
    -u +15 -U +15 \
    -a AGATCGGAAGAGCACAC -A AGATCGGAAGAGCGTCG \
    --minimum-length=30 \
    --output=1-trim/${stem}-R1.fastq.gz \
    --paired-output=1-trim/${stem}-R2.fastq.gz \
    $input1 $input2
