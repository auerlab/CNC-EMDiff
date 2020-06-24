#!/bin/sh -e

##########################################################################
#   Script description:
#       Align reads to reference genome
#
#   Usage:
#       SLURM cluster:
#           sbatch 4-bwa-mem.sbatch
#       No cluster:
#           [env SLURM_CPUS_PER_TASK=n] \
#               ./slurm-sim 4-bwa-mem.sbatch 1 18 |& tee 4.log
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

# SLURM parameters
# Not yet sure how many threads are useful.  Maybe use --exclusive and all
# available cores? Smaller thread counts are usually more efficient as
# they will require less communication overhead.  If cores are limited,
# this will get the job done faster.
#SBATCH --array=1-18 --cpus-per-task=2
# Based on 2 threads on Albacore cluster
# 12 threads on Peregrine used close to 6g
#SBATCH --mem=5g
#SBATCH --output=4-bwa-mem/slurm-%A_%a.out
#SBATCH --error=4-bwa-mem/slurm-%A_%a.err

if [ -z $SLURM_ARRAY_TASK_ID ]; then
    printf "$0 must be run under sbatch or ./slurm-sim\n"
    exit 1
fi

# Defaults if not set by SLURM or user
: ${SLURM_CPUS_PER_TASK:=2}

# Document version of software used for publication
bwa > 4-bwa-mem/bwa-version.txt 2>&1

# One iteration if running under SLURM, all iterations otherwise
cd 4-bwa-mem
printf "Running task $SLURM_ARRAY_TASK_ID...\n"

infile1=../1-trim/*_S${SLURM_ARRAY_TASK_ID}_L001-R1.fastq.gz
infile2=../1-trim/*_S${SLURM_ARRAY_TASK_ID}_L001-R2.fastq.gz
outfile=$(echo $infile1 | cut -d / -f 3 | sed -e 's|-R1.fastq.gz|.sam|')
ref_file=$(../Reference/ref-file-name.sh)

srun bwa mem -M -t $SLURM_CPUS_PER_TASK \
    ../Reference/$ref_file $infile1 $infile2 > $outfile
