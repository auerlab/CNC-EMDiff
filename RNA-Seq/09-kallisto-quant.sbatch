#!/bin/sh -e

##########################################################################
#   Description:
#       Run kallisto quantification for each RNA sample.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#       
#   History:
#   Date        Name        Modification
#   2019-11-??  Jason Bacon Begin
##########################################################################

#SBATCH --array=1-15
#SBATCH --cpus-per-task=2
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem=6000
#SBATCH --output=Logs/09-kallisto-quant/slurm-%A_%a.out
#SBATCH --error=Logs/09-kallisto-quant/slurm-%A_%a.err

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

# Document software versions used for publication
uname -a
kallisto version
pwd

gtf=$(Reference/gtf-filename.sh)

# If using hdf5, you may need this:
# https://github.com/pachterlab/kallisto/issues/197
# export HDF5_USE_FILE_LOCKING=FALSE

# 6-merge-bams.sbatch relies on sample N being in Data/09-kallisto-quant/N
# The sample number comes after -sample in the filename, e.g.
# chondro-sample4-rep2-time1-R1.fastq.xz is sample 4

# kallisto 0.46.1 can't handle xz and will simply seg fault rather than
# issue an error message.  If your trimmed fastq files are in xz format,
# this will convert to gzip format.
# Convert xz to gz rather than raw to reduce NFS load from compute nodes
xz1=$(echo Data/04-trim/*sample${SLURM_ARRAY_TASK_ID}-*-R1.fastq.xz)
xz2=$(echo Data/04-trim/*sample${SLURM_ARRAY_TASK_ID}-*-R2.fastq.xz)
gz1=$(echo Data/04-trim/*sample${SLURM_ARRAY_TASK_ID}-*-R1.fastq.gz)
gz2=$(echo Data/04-trim/*sample${SLURM_ARRAY_TASK_ID}-*-R2.fastq.gz)
test -e $gz1 || xzcat $xz1 | gzip --stdout > $gz1
test -e $gz2 || xzcat $xz2 | gzip --stdout > $gz2

# Kallisto requires an output subdirectory for each sample
stem=$(basename ${gz1%-R1*})
out_dir=Data/09-kallisto-quant/$stem
mkdir -p $out_dir

# Manual says a GTF is needed.  Kallisto aborts using GFF3.
set -x
kallisto quant \
    --genomebam \
	--gtf=Data/07-reference/$gtf \
	--chromosomes=Data/07-reference/chromosome-sizes.tsv \
    --bootstrap-samples=500 \
    --threads=2 \
    --index=Data/08-kallisto-index/all-but-xy.index \
    --output-dir=$out_dir $gz1 $gz2

