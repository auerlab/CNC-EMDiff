#!/bin/sh -e

#SBATCH --mem=3g

ref_file=$(./ref-file-name)
if [ ! -e $ref_file ]; then
    cd Reference
    ./fetch-ref.sh
    ./concat-ref.sh
    cd ..
fi
bwa index $ref_file
