#!/bin/sh -e

# Don't clobber previous results
for dir in 1-trim 2-qc 4-kallisto-quant; do
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
trim_job_id=$(sbatch 1-trim.sbatch | awk '{ print $4 }')

# Run fastqc quality checks on raw and trimmed samples
sbatch --dependency=afterok:$trim_job_id 2-qc.sbatch

# Run kallisto index if necessary. This is independent of the samples, so
# it may be left alone between quantification runs.
if [ ! -e 3-kallisto-index/all-but-xy.index ]; then
    index_job_id=$(sbatch 3-kallisto-index.sbatch | awk '{ print $4 }')
    quant_dependency="--dependency=afterok:${trim_job_id}:$index_job_id"
else
    quant_dependency="--dependency=afterok:$trim_job_id"
fi

# Kallisto quantification.  Requires trimming and indexing to be completed.
quant_job_id=$(sbatch $quant_dependency 4-kallisto-quant.sbatch \
    | awk '{ print $4 }')

sbatch --dependency=afterok:$quant_job_id 5-merge-bams.sbatch
