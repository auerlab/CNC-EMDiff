#!/bin/sh -e

date=$(date +%Y-%m-%d-%H:%M)

# Inefficient, but portable way to get last argument (filename)
for last_arg in "$@"; do
    echo '' > /dev/null
done
filename="$last_arg"

script=$1
shift

base=$(basename $script)
log_dir=Logs/${base%.sh}
base=$(basename $filename)
output_log=$log_dir/$base-$date.out
error_log=$log_dir/$base-$date.err

cmd="$script $@ > $output_log 2> $error_log"
printf "===> $cmd\n"
eval $cmd
