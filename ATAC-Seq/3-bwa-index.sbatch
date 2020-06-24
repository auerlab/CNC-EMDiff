#!/bin/sh -e

##########################################################################
#   Script description:
#       Generate index for BWA aligner
#
#   Usage:
#       sbatch 3-bwa-index.sbatch
#       ./3-bwa-index.sbatch |& tee 3.log
#
#   History:
#       Based on the work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#   Date        Name        Modification
#   2020-02-10  Jason Bacon Begin
##########################################################################

##########################################################################
#   Main
##########################################################################

# SLURM parameters
#SBATCH --mem=3g
#SBATCH --output=3-bwa-index/slurm-%A_%a.out
#SBATCH --error=3-bwa-index/slurm-%A_%a.err

if [ -z $SLURM_ARRAY_TASK_ID ]; then
    printf "$0 must be run under sbatch or ./slurm-sim\n"
    exit 1
fi

bwa > Reference/bwa-version.txt 2>&1

cd Reference
ref_file=$(./ref-file-name.sh)
if [ ! -e $ref_file ]; then
    ./fetch-ref.sh
    ./concat-ref.sh
    bwa index $ref_file
fi
