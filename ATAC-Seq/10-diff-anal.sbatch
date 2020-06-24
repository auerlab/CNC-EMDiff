#!/bin/sh -e

##########################################################################
#   Script description:
#       Find differential accessibility across time points
#
#   Usage:
#       SLURM Cluster:
#           sbatch 10-diff-anal.sbatch
#       No cluster:
#           ./slurm-sim 10-diff-anal.sbatch 1 2 |& tee 1.log
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

# SLURM parameters
#SBATCH --array=1-2
#SBATCH --mem=3g
#SBATCH --output=10-diff-anal/slurm-%A_%a.out
#SBATCH --error=10-diff-anal/slurm-%A_%a.err

if [ -z $SLURM_ARRAY_TASK_ID ]; then
    printf "$0 must be run under sbatch or ./slurm-sim\n"
    exit 1
fi

case $SLURM_ARRAY_TASK_ID in
1)
    cell_type=CCA
    ;;
2)
    cell_type=NCA
    ;;
esac

# Document software versions used for publication
uname -a > 10-diff-anal/os-version.txt 2>&1
R --version > 10-diff-anal/R-version.txt 2>&1

cd 10-diff-anal
../10-diff-anal.R $cell_type
