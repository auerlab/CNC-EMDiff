#!/bin/sh -e

##########################################################################
#   Script description:
#       Generate index for BWA aligner
#
#   Usage:
#       sbatch 07-bwa-index.sbatch
#       ./07-bwa-index.sbatch |& tee 3.log
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

#SBATCH --mem=3g
#SBATCH --output=Logs/07-bwa-index/slurm-%A_%a.out
#SBATCH --error=Logs/07-bwa-index/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
bwa || true
pwd

genome_file=$(../RNA-Seq/Reference/genome-filename.sh)

genome_dir=../RNA-Seq/Data/07-reference
if [ ! -e $genome_dir/$genome_file ]; then
    save_cwd=$(pwd)
    cd ../RNA-Seq
    Reference/build-genome.sh
    cd $save_cwd
fi

cd Data/07-bwa-index
ln -s ../../$genome_dir/$genome_file
bwa index $genome_file
