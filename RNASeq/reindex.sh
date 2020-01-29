#!/bin/sh -e

cat << EOM

This will remove the contents of 3-kallisto-index and rebuild them using
3-kallisto-index.sbatch.

EOM
read -p "Are you sure you want to do this? yes/[no] " sure
if [ 0"$sure" = 0yes ]; then
    rm -rf 3-kallisto-index/*
    sbatch 3-kallisto-index.sbatch
    squeue
fi
