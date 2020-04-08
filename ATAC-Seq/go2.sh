#!/bin/sh -e

# Don't clobber previous results
for dir in 3-bwa-index 4-bwa-mem 5-qc-sam 6-remove-duplicates 7-macs-peaks \
	   8-merge-bams; do
    if [ -e $dir ]; then
	pre_existing="$pre_existing $dir"
    fi
done
if [ 0"$pre_existing" != 0 ]; then
    printf "The following directories already exist: $pre_existing\n"
    printf "Remove or rename them and try again.\n"
    exit 1
fi

if [ ! -e 2-qc/Trimmed/NCA3C_S18_L001-R2_fastqc.zip ]; then
    cat << EOM

You must perform trimming and manually examine quality checks first.
Run go1.sh or manually run 1-trim.sbatch and 2-qc.sbatch.

EOM
    exit 1
fi

# Some dirs referenced by #SBATCH directives, so cannot be created by
# sbatch scripts
./0-mkdirs

# Generate genome index for BWA
bwa_index_job_id=$(sbatch "$@" 3-bwa-index.sbatch | awk '{ print $4 }')

# Align reads to genome
bwa_mem_job_id=$(sbatch "$@" --dependency=afterok:$bwa_index_job_id \
		 4-bwa-mem.sbatch | awk '{ print $4 }')

# QC alignments
sbatch "$@" --dependency=afterok:$bwa_mem_job_id 5-qc-sam.sbatch

remove_dups_job_id=$(sbatch "$@" --dependency=afterok:$bwa_mem_job_id \
		     6-remove-duplicates.sbatch | awk '{ print $4 }')

macs_job_id=$(sbatch "$@" --dependency=afterok:$remove_dups_job_id \
		     7-macs-peaklets.sbatch | awk '{ print $4 }')

merge_bams_job_id=$(sbatch "$@" --dependency=afterok:$remove_dups_job_id \
		     8-merge-bams.sbatch | awk '{ print $4 }')

process_peaks_job_id=$(sbatch "$@" --dependency=afterok:$macs_job_id \
		       9-process-peaks.sbatch | awk '{ print $4 }')

# Create an sbatch script to run the R script on both cell types
sbatch "$@" --dependency=afterok:$process_peaks_job_id 10-diff-anal.sbatch
