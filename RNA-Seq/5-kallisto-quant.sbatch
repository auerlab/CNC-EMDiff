#!/bin/sh -e

#SBATCH --mem-per-cpu=4000
#SBATCH --array=1-18
#SBATCH --output=Data/4-kallisto-quant/slurm-%A_%a.out
#SBATCH --error=Data/4-kallisto-quant/slurm-%A_%a.err
##SBATCH --exclude=compute-[001-008,012]
##SBATCH --partition=batch,256g

# Run kallisto with 500 bootstraps for Sleuth
# --genomebam is needed to generate a genome-mapped BAM file for browsing with
# IGV.  It requires --gtf and --chromosomes. --chromosomes requires a TSV file
# with chromosome name and length on each line.  The name must match the names
# used in the GTF.
# https://github.com/pachterlab/kallisto/issues/155
# The format and source of the chromosomes TSV is not clear.  I generated one
# using an Ensemble GFF with Reference/create-chrom-sizes.sh

mkdir -p Data/4-kallisto-quant/$SLURM_ARRAY_TASK_ID
set -x
kallisto quant \
    --genomebam \
	--gtf=Reference/Mus_musculus.GRCm38.98.gtf \
	--chromosomes=Reference/chromosome-sizes.tsv \
    -i Data/3-kallisto-index/all-but-xy.index \
    -o Data/4-kallisto-quant/$SLURM_ARRAY_TASK_ID \
    -b 500 \
    1-trim/*_S${SLURM_ARRAY_TASK_ID}_L002-R1.fastq.gz \
    1-trim/*_S${SLURM_ARRAY_TASK_ID}_L002-R2.fastq.gz
