#!/bin/sh -e

#SBATCH --mem=3g
#SBATCH --output=3-bwa-index/slurm-%A_%a.out
#SBATCH --error=3-bwa-index/slurm-%A_%a.err

cd Reference
bwa > bwa-version.txt 2>&1
ref_file=$(./ref-file-name.sh)
if [ ! -e $ref_file ]; then
    ./fetch-ref.sh
    ./concat-ref.sh
    bwa index $ref_file
fi
