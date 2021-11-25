#!/bin/sh -e

#SBATCH --mem=6000
#SBATCH --array=1-18
#SBATCH --ntasks=2
#SBATCH --output=Logs/5-kallisto-quant/slurm-%A_%a.out
#SBATCH --error=Logs/5-kallisto-quant/slurm-%A_%a.err

##############################################################################
# Run kallisto with 500 bootstraps for Sleuth
#
# --genomebam is needed to generate a genome-mapped BAM file for browsing with
# IGV.  It requires --gtf and --chromosomes. --chromosomes requires a TSV file
# with chromosome name and length on each line.  The chromosome names in the
# TSV must exactly match the names in the GTF.
# https://github.com/pachterlab/kallisto/issues/155
#
# The format and source of the chromosomes TSV is not clearly documented.
# I generated one using an Ensemble GFF with Reference/create-chrom-sizes.sh.
# GTF does not contain chromosome features.

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}

# Document software versions used for publication
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
    uname -a > Logs/5-kallisto-quant/os-version-$SLURM_JOB_ID.txt 2>&1
    kallisto version > Logs/5-kallisto-quant/kallisto-version-$SLURM_JOB_ID.txt 2>&1
fi

build=$(../Common/genome-build.sh)
release=$(../Common/genome-release.sh)

# If using hdf5, you may need this:
# https://github.com/pachterlab/kallisto/issues/197
# export HDF5_USE_FILE_LOCKING=FALSE

# 6-merge-bams.sbatch relies on sample N being in Data/5-kallisto-quant/N
# The sample number comes after -sample in the filename, e.g.
# chondro-sample4-rep2-time1-R1.fastq.xz is sample 4

# kallisto 0.46.1 can't handle xz and will simply seg fault rather than
# issue an error message.
# Convert xz to gz rather than raw to reduce NFS load from compute nodes
xz1=$(echo Data/1-trim/*-sample${SLURM_ARRAY_TASK_ID}-*-R1.fastq.xz)
xz2=$(echo Data/1-trim/*-sample${SLURM_ARRAY_TASK_ID}-*-R2.fastq.xz)
gz1=${xz1%.xz}.gz
gz2=${xz2%.xz}.gz
test -e $gz1 || xzcat $xz1 | gzip --stdout > $gz1
test -e $gz2 || xzcat $xz2 | gzip --stdout > $gz2

# Kallisto requires an output subdirectory for each sample
stem=$(basename ${xz1%-R1*})
out_dir=Data/5-kallisto-quant/$stem
echo $out_dir
mkdir -p $out_dir

# Manual says a GTF is needed.  Kallisto aborts using GFF3.
set -x
kallisto quant \
    --genomebam \
	--gtf=Data/3-reference/Mus_musculus.GRCm$build.$release.chr.gtf \
	--chromosomes=Data/3-reference/chromosome-sizes.tsv \
    --bootstrap-samples=500 \
    --threads=2 \
    --index=Data/4-kallisto-index/all-but-xy.index \
    --output-dir=$out_dir $gz1 $gz2

# Remove temporary .gz files, keep .xz which is ~40% smaller
rm -f $gz1 $gz2
