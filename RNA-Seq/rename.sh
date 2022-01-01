#!/bin/sh -e

usage()
{
    printf "Usage: $0 high low increment\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 3 ]; then
    usage
fi
high=$1
low=$2
inc=$3

c=$high
while [ $c -ge $low ]; do
    prefix=$(printf "%02d" $c)
    script=$(ls $prefix-*)
    base=${script%.*}
    new_prefix=`printf "%02d" $((c + inc))`
    new_script=$(echo $script | sed -e "s|$prefix|$new_prefix|")
    new_base=$(echo $base | sed -e "s|$prefix|$new_prefix|")
    printf "$base -> $new_base\n"
    
    sed -i '' -e "s|$base|$new_base|g" *.sbatch *.sh *.R
    git mv $script $new_script
    mv Data/$base Data/$new_base || true
    mv Logs/$base Logs/$new_base || true
    
    c=$((c - 1))
done
