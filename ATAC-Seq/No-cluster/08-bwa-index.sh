#!/bin/sh -e

##########################################################################
#   Script description:
#       Generate index for BWA aligner
#
#   Usage:
#       sbatch 08-bwa-index.sbatch
#       ./08-bwa-index.sbatch |& tee 3.log
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

# Document software versions used for publication
uname -a
bwa || true
pwd

genome_file=$(../../RNA-Seq/Reference/genome-filename.sh)

genome_dir=../../RNA-Seq/Results/07-reference
if [ ! -e $genome_dir/$genome_file ]; then
    save_cwd=$(pwd)
    cd ../../RNA-Seq
    Reference/build-genome.sh
    cd $save_cwd
fi

set -x
cd Results/08-bwa-index
ln -sf ../../$genome_dir/$genome_file
bwa index $genome_file
