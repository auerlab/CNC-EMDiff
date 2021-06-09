#!/bin/sh -e

cat << EOM

This will remove the contents of 6-merge-bams and rebuild them using
3-build-reference.sbatch.

EOM
read -p "Are you sure you want to do this? yes/[no] " sure
if [ 0"$sure" = 0yes ]; then
    rm -rf Data/6-merge-bams/*
    sbatch 6-merge-bams.sbatch
    squeue
fi
