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

# Document software versions used for publication
uname -a > 3-bwa-index/os-version.txt 2>&1
bwa 2>&1 | fgrep Version > 3-bwa-index/bwa-version.txt || true

cd 3-bwa-index
ref_file=$(../ref-file-name.sh)
if [ ! -e $ref_file ]; then
    ../fetch-ref.sh
    ../concat-ref.sh
    bwa index $ref_file
fi
