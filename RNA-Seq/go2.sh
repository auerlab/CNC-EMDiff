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
for dir in 4-kallisto-quant; do
    if [ -e $dir ]; then
	pre_existing="$pre_existing $dir"
    fi
done
if [ 0"$pre_existing" != 0 ]; then
    printf "The following directories already exist: $pre_existing\n"
    printf "Remove or rename them and try again.\n"
    exit 1
fi

# FIXME: Check for output files
if [ ! -e 2-qc/Trimmed ]; then
    cat << EOM

You must trim and check quality before running this script.

EOM
    exit 1
fi

# Some dirs referenced by #SBATCH directives, so cannot be created by
# sbatch scripts
./0-mkdirs

# Run kallisto index if necessary. This is independent of the samples, so
# it may be left alone between quantification runs.
if [ ! -e 3-kallisto-index/all-but-xy.index ]; then
    index_job_id=$(sbatch "$@" 3-kallisto-index.sbatch | awk '{ print $4 }')
    quant_dependency="--dependency=afterok:$index_job_id"
fi

# Kallisto quantification.  Requires trimming and indexing to be completed.
quant_job_id=$(sbatch "$@" $quant_dependency 4-kallisto-quant.sbatch \
    | awk '{ print $4 }')

sbatch "$@" --dependency=afterok:$quant_job_id 5-merge-bams.sbatch
