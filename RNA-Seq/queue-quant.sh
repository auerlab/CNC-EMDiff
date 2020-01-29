#!/bin/sh -e

##########################################################################
#   Script description:
#       Queue kallisto quantification for after trim job
#       
#   History:
#   Date        Name        Modification
#   2019-10-09  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 trim-job-id\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
jobid=$1

cat << EOM

This will remove the contents of 4-kallisto-quant and rebuild them using
4-kallisto-quant.sbatch.

EOM
read -p "Are you sure you want to do this? yes/[no] " sure
if [ 0"$sure" = 0yes ]; then
    rm -rf 4-kallisto-quant/*
    sbatch --dependency=afterok:$jobid 4-kallisto-quant.sbatch
    squeue
fi
