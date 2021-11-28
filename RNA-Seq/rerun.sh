#!/bin/sh -e

##########################################################################
#   Description:
#       Remove output files and logs from a previous run and resubmit
#       
#   History:
#   Date        Name        Modification
#   2021-11-27  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 script.sbatch\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
script=$1

if ! echo $script | grep -q "\.sbatch$"; then
    printf "$0 is only for .sbatch scripts.\n"
    usage
fi

base=${script%.sbatch}
printf "Are you sure you want to remove results from Data/$base? y/[n] "
read sure
if [ 0"$sure" = 0y ]; then
    rm -r Data/$base
    ./00-organize.sh
    sbatch $script
fi
