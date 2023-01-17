#!/bin/sh -e

for script in [012]*-*.s*; do
    dir=Logs/${script%.*}
    if [ -d $dir ]; then
	printf "%-40s: " $script
	if [ -z "$(ls $dir)" ]; then
	    printf "Not submitted.\n"
	else
	    printf "Submitted\n"
	fi
    fi
done
