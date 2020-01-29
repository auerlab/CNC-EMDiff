#!/bin/sh -e

cat << EOM

This will remove the contents of 4-kallisto-quant and rebuild them using
4-kallisto-quant.sbatch.

EOM
read -p "Are you sure you want to do this? yes/[no] " sure
if [ 0"$sure" = 0yes ]; then
    rm -rf 4-kallisto-quant/*
    sbatch 4-kallisto-quant.sbatch
    squeue
fi
