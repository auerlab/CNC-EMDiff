#!/bin/sh -e

##########################################################################
#   Script description:
#       Consolidate FastQC reports into one for easy viewing
#
#   Dependencies:
#       Requires raw FastQC results.  Run after *-qc-raw.sh.
##########################################################################

##########################################################################
# multiqc: LC_ALL and LANG must be set to a UTF-8 character set
# in your environment in order for the click module to function.
# export LC_ALL=en_US.UTF8

LC_ALL=en_US.UTF-8
LANG=en_US.utf-8
export LC_ALL LANG

date=$(date +%Y-%m-%d-%H:%M)
set -x
multiqc --outdir Results/04-multiqc-raw Results/03-qc-raw \
    2>&1 | tee Logs/04-multiqc-raw/$date.out
