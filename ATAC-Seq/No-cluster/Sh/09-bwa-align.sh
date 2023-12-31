#!/bin/sh -e

##########################################################################
#   Script description:
#       Align reads to reference genome
#
#   Usage:
#       SLURM cluster:
#           sbatch 09-bwa-mem.sbatch
#
#   History:
#       Based on the work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#   Date        Name        Modification
#   2020-02-20  Jason Bacon Begin
##########################################################################

##########################################################################
#   Main
##########################################################################

if [ $# != 2 ]; then
    printf "Usage: $0 threads-per-job input fasta-R1\n"
    exit 1
fi
threads_per_job=$1
infile1=$2
infile2=$(echo $infile1 | sed -e 's|R1|R2|')

# Document software versions used for publication
uname -a
bwa || true
pwd

ref_file=$(../../RNA-Seq/Reference/genome-filename.sh)

# One iteration if running under SLURM, all iterations otherwise
cd Results/09-bwa-align
outfile=$(basename ${infile1%-R1.fastq.gz}).sam

set -x
pipe1=$(basename $infile1 | sed -e 's|zst|fifo|')
pipe2=$(basename $infile2 | sed -e 's|zst|fifo|')
rm -f $pipe1 $pipe2
mkfifo $pipe1 $pipe2
zstdcat ../../$infile1 > $pipe1 &
zstdcat ../../$infile2 > $pipe2 &

date
set -x
bwa mem -M -t $threads_per_job \
    ../08-bwa-index/$ref_file $pipe1 $pipe2 > $outfile
set +x
wait
date
rm -f $pipe1 $pipe2
