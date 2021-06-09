#!/bin/sh -e

cat << EOM

This will remove the contents of 2-qc and rebuild them using 2-qc.sbatch.

EOM
read -p "Are you sure you want to do this? yes/[no] " sure
if [ 0"$sure" = 0yes ]; then
    rm -rf Data/2-qc/*
    sbatch 2-qc.sbatch
    squeue
fi
