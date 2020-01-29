#!/bin/sh -e

##########################################################################
#   Script description:
#       Simulate SLURM env by setting SLURM_ARRAY_TASK_ID and running
#       sbatch scripts in a loop.  This allows using the sbatch scripts
#       outside the SLURM env, such as on a development server or
#       workstation.
#
#       It is advisable to run this script under nohup, screen, or other
#       tool that allows it to continue after a dropped connection.
#
#       E.g. nohup ./slurm-sim.sh 1-trim >& slurm-sim.out &
#
#   Arguments:
#       Analysis step name (should be same as stem name of sbatch script)
#       
#   History:
#   Date        Name        Modification
#   2019-09-29  Jason Wayne BaconBegin
##########################################################################

usage()
{
    printf "Usage: $0 step-name (e.g. 1-trim or 2-qc)\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi

step_name="$1"
if [ ! -e $step_name.sbatch ]; then
    printf "Invalid step-name: $step_name (no $step_name.sbatch found)\n"
    exit 1
fi

mkdir -p $step_name
c=1
set -x
while [ $c -le 18 ]; do
    export SLURM_ARRAY_TASK_ID=$c
    sh "$step_name.sbatch" | tee $step_name/$c.out 2| tee $step_name/$c.err
    c=$((c + 1))
done
