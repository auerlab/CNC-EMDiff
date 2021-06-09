#!/usr/bin/env bash

#SBATCH --mem=4000
#SBATCH --array=1-18
#SBATCH --output=Data/5-kallisto-quant/slurm-%A_%a.out
#SBATCH --error=Data/5-kallisto-quant/slurm-%A_%a.err

set -e

##############################################################################
# Run kallisto with 500 bootstraps for Sleuth
# --genomebam is needed to generate a genome-mapped BAM file for browsing with
# IGV.  It requires --gtf and --chromosomes. --chromosomes requires a TSV file
# with chromosome name and length on each line.  The name must match the names
# used in the GTF.
# https://github.com/pachterlab/kallisto/issues/155
# The format and source of the chromosomes TSV is not clear.  I generated one
# using an Ensemble GFF with Reference/create-chrom-sizes.sh

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}

# Document software versions used for publication
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
    uname -a > Data/5-kallisto-quant/os-version-$SLURM_JOB_ID.txt 2>&1
    kallisto version > Data/5-kallisto-quant/kallisto-version-$SLURM_JOB_ID.txt 2>&1
fi

build=$(../Common/genome-build.sh)
release=$(../Common/genome-release.sh)

# https://github.com/pachterlab/kallisto/issues/197
export HDF5_USE_FILE_LOCKING=FALSE

# Need bash for process substitution: <()
kallisto quant \
    --genomebam \
	--gtf=Data/3-reference/Mus_musculus.GRCm$build.$release.gtf \
	--chromosomes=Data/3-reference/chromosome-sizes.tsv \
    -i Data/4-kallisto-index/all-but-xy.index \
    -o Data/5-kallisto-quant \
    -b 500 \
    <( xzcat Data/1-trim/*_S${SLURM_ARRAY_TASK_ID}_L002-R1.fastq.xz ) \
    <( xzcat Data/1-trim/*_S${SLURM_ARRAY_TASK_ID}_L002-R2.fastq.xz )
