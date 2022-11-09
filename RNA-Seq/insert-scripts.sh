#!/bin/sh -e

usage()
{
    printf "Usage: $0 start-index\n"
    exit 1
}


##########################################################################
#   Function description:
#       
#   Arguments:
#       
#   Returns:
#       
#   History:
#   Date        Name        Modification
#   2022-11-09  Jason Bacon Begin
##########################################################################

bump_prefix()
{
    if [ $# != 1 ]; then
	printf "Usage: bump_prefix filename\n"
	exit 1
    fi
    original_script=$1
    #echo $original_script
    original_prefix=${original_script%%-*}
    #echo $original_prefix
    new_prefix=`printf "%02d" $(($original_prefix + 1))`
    new_script=$new_prefix-${original_script#*-}
    original_base=${original_script%.*}
    new_base=${new_script%.*}
    printf "$original_base -> $new_base\n"
    
    sed -i '' -e "s|$original_base|$new_base|g" *.sbatch *.sh *.R
    git mv $original_script $new_script
    if [ -e Data/$original_base ]; then
	mv Data/$original_base Data/$new_base
    fi
    if [ -e Logs/$original_base ]; then
	mv Logs/$original_base Logs/$new_base || true
    fi
    printf "\n"
    return 0
}

##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
start=$1

c=$start
while [ -e $c-*.sbatch ] || [ -e $c-*.sh ]; do
    #ls $c-*.sbatch
    c=$(($c + 1))
done

while [ $c -gt $start ]; do
    original_index=$(($c - 1))
    prefix=$(printf "%02d" $original_index)
    
    if [ -e $original_index-*.sh ]; then
	bump_prefix $original_index-*.sh
    elif [ -e $original_index-*.sbatch ]; then
	bump_prefix $original_index-*.sbatch
    fi
    
    c=$((c - 1))
done
