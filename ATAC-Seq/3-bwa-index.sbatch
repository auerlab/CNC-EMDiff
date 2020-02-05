#!/bin/sh -e

#SBATCH --mem=3g

ref_file=$(./ref-file-name)

bwa index $ref_file
