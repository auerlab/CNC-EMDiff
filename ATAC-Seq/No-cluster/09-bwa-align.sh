#!/bin/sh -e

##########################################################################
#   Script description:
#       Use hisat2 to align reads to a genome reference.
#
#   Dependencies:
#       Requires trimmed reads and a reference genome.  Run after
#       *-trim.sh and *-reference.sh.
##########################################################################

set -x
hw_threads=$(../../Common/get-hw-threads.sh)
hw_mem=$(../../Common/get-hw-mem.sh)
hw_gib=$(( $hw_mem / 1024 / 1024 / 1024 ))

# For xenopus, hisat2 jobs take about 4.3 GB
# For axolotl, quite a bit more
if pwd | fgrep XenONReg; then
    jobs=$(( $hw_gib / 6 ))
else
    jobs=$(( $hw_gib / 16 ))
fi

# Not enough RAM according to calculations above, but try anyway
if [ $jobs = 0 ]; then
    jobs=1
fi
threads_per_job=$(( $hw_threads / $jobs ))

# Use at least 4 thread per job, since the jobs will run almost 4 times
# as fast, and fewer jobs means less disk contention.  If CPU utilization
# is less than 90% per core (360% for a 4-thread job), you probably have
# too much disk contention and need to reduce the number of jobs.
if [ $threads_per_job -lt 4 ]; then
    threads_per_job=4
    jobs=$(( $hw_threads / $threads_per_job ))
fi

# Tried GNU parallel and ran into bugs.  Xargs just works.
ls Results/05-trim/*-R1.fastq.zst | \
    xargs -n 1 -P $jobs ../../Common/redirect.sh \
    Xargs/09-bwa-align.sh $threads_per_job
