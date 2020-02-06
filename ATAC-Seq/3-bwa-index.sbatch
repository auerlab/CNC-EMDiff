#!/bin/sh -e

#SBATCH --mem=3g
#SBATCH --output=Reference/slurm-%A_%a.out
#SBATCH --error=Reference/slurm-%A_%a.err

ref_file=$(./ref-file-name)
if [ ! -e $ref_file ]; then
    cd Reference
    ./fetch-ref.sh
    ./concat-ref.sh
    cd ..
fi
bwa index $ref_file
