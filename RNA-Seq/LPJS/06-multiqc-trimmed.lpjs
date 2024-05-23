#!/bin/sh -e

#########################################################################
#   MultiQC is optional, but helpful for visualizing raw read quality
#   and trimming results
#
#   Dependencies:
#       Requires raw FastQC results.  Run after *-qc-raw.sbatch.

#lpjs jobs 1
#lpjs procs-per-job 1
#lpjs min-procs-per-node procs-per-job
#lpjs pmem-per-proc 100MiB
#lpjs log-dir Logs/03-multiqc-raw

# multiqc: LC_ALL and LANG must be set to a UTF-8 character set
# in your environment in order for the click module to function.
export LC_ALL=en_US.UTF-8

cd Results/06-multiqc-trimmed
rm -rf *
multiqc --version > ../../Logs/06-multiqc-trimmed/multiqc-version.txt 2>&1
multiqc ../05-qc-trimmed 2>&1 | tee ../../Logs/06-multiqc-trimmed/multiqc.out
