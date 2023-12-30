#!/bin/sh -e

##########################################################################
#   Script description:
#       Trim adapters, poly-A tails, and low-quality bases from reads.
#
#   Dependencies:
#       Requires directory structure.  Run after *-organize.sh.
##########################################################################

hw_threads=$(../../Common/get-hw-threads.sh)
jobs=$(($hw_threads / 2))

# Limit jobs to avoid overwhelming disk I/O.
# Reduce this if using a single disk, increase it if using a fast
# RAID or SSD.  The "top" command should show fastq-trim processes
# using 80% or more CPU each.  If it's lower than this, you probably
# have too many jobs competing for disk.
max_jobs_for_disk=5
if [ $jobs -gt $max_jobs_for_disk ]; then
    jobs=$max_jobs_for_disk
fi

# Tried GNU parallel and ran into bugs.  Xargs just works.
ls Results/01-organize/Raw-renamed/*-R1.fastq.xz | \
    xargs -n 1 -P $jobs ../../Common/redirect.sh Sh/05-trim.sh
