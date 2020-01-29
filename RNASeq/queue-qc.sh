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

This will remove the contents of 2-qc and rebuild them using
2-qc.sbatch.

EOM
read -p "Are you sure you want to do this? yes/[no] " sure
if [ 0"$sure" = 0yes ]; then
    rm -rf 2-qc/*
    sbatch --dependency=afterok:$jobid 2-qc.sbatch
    squeue
fi
