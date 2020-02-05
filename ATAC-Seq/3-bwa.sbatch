#!/bin/sh -e

#SBATCH --array=1-18
# Really?
#SBATCH --mem-per-cpu=64000

sample=( 0ATAC-1 0ATAC-2 0ATAC-3 12ATAC-1 12ATAC-2 12ATAC-3 2ATAC-1 2ATAC-2 2ATAC-3 4ATAC-1 4ATAC-2 4ATAC-3 7ATAC-2 7ATAC-3 )
bwa mem -M ~/Data/Optic-Regen/tools/REF_DATA/danio_trans.fa \
    ${homedir}TRIMMED/${s}_R1_val_1.fq.gz \
    ${homedir}TRIMMED/${s}_R2_val_2.fq.gz \
    > ${homedir}ALIGNED/${s}_trans.sam"
done
