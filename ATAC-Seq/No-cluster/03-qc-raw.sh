#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on raw data.
#
#   Dependencies:
#       Requires directory structure.  Run after *-organize.sh.
##########################################################################

# FastQC can utilize 2 cores, so divide total hardware cores by 2
# to determine the number of simultaneous jobs
hw_threads=$(../../Common/get-hw-threads.sh)
jobs=$(($hw_threads / 2))

# Tried GNU parallel and ran into bugs.  Xargs just works.
date
ls Results/01-organize/Raw-renamed/*.fastq.xz | \
    xargs -n 1 -P $jobs ../../Common/redirect.sh Sh/03-qc-raw.sh
date
