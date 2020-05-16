#!/bin/sh -e

##########################################################################
#   Script description:
#       Simulate SLURM array environment by running an sbatch script
#       in a loop with SLURM_ARRAY_TASK_ID set for each iteration.
#       
#   History:
#   Date        Name        Modification
#   2020-05-16  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 script first-index last-index\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 3 ]; then
    usage
fi

script=$1
first=$2
last=$3

SLURM_ARRAY_TASK_ID=$first
while [ $SLURM_ARRAY_TASK_ID -le $last ]; do
    export SLURM_ARRAY_TASK_ID
    ./$script
done
