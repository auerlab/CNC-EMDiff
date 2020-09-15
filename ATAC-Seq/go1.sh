#!/bin/sh -e

cat << EOM

This will automatically run multiple stages of the pipeline without meaningful
validation of intermediate outputs.

This should only be done if the latest pipeline changes have been well-tested
and you are certain that everything is working.

EOM
read -p "Continue?  yes/[no] " continue
if [ 0$continue != 0yes ]; then
    printf "OK, aborting.\n"
    exit
fi

# Don't clobber previous results
for dir in 1-trim 2-qc; do
    if [ -e $dir ]; then
	pre_existing="$pre_existing $dir"
    fi
done
if [ 0"$pre_existing" != 0 ]; then
    printf "The following directories already exist: $pre_existing\n"
    printf "Remove or rename them and try again.\n"
    exit 1
fi

# Some dirs referenced by #SBATCH directives, so cannot be created by
# sbatch scripts
./0-mkdirs

# Trim adapters and low-quality data
trim_job_id=$(sbatch "$@" 1-trim.sbatch | awk '{ print $4 }')

sbatch "$@" --dependency=afterok:$trim_job_id 2-qc.sbatch
