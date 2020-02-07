#!/bin/sh -e

#SBATCH --mem=3g
#SBATCH --output=Reference/slurm-%A_%a.out
#SBATCH --error=Reference/slurm-%A_%a.err

cd Reference
ref_file=$(./ref-file-name.sh)
if [ ! -e $ref_file ]; then
    ./fetch-ref.sh
    ./concat-ref.sh
    bwa index $ref_file
fi

