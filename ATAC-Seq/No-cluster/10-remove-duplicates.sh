#!/bin/sh -e

##########################################################################
#   Script description:
#       Use BWA to align reads to a genome reference.
#
#   Dependencies:
#       Requires trimmed reads and a reference genome.  Run after
#       *-trim.sh and *-reference.sh.
##########################################################################

set -x
hw_threads=$(../../Common/get-hw-threads.sh)
hw_mem=$(../../Common/get-hw-mem.sh)
hw_gib=$(( $hw_mem / 1024 / 1024 / 1024 ))

# Use lesser of cores and mem / 2 GiB
# FIXME: Really only needs about 2.  This is for testing to reserve
# mem for other tasks.
jobs=$(( $hw_gib / 8 ))
if [ $hw_threads -lt $jobs ]; then
    jobs=$hw_threads
fi
threads_per_job=1
printf "jobs = $jobs\n"

# Tried GNU parallel and ran into bugs.  Xargs just works.
ls Results/09-bwa-align |
    xargs -n 1 -P $jobs ../../Common/redirect.sh \
    Xargs/10-remove-duplicates.sh
