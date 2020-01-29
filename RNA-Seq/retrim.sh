#!/bin/sh -e

cat << EOM

This will remove the contents of 1-trim and rebuild them using 1-trim.sbatch.

EOM
read -p "Are you sure you want to do this? yes/[no] " sure
if [ 0"$sure" = 0yes ]; then
    rm -rf 1-trim/*
    sbatch 1-trim.sbatch
    squeue
fi
