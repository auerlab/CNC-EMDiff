#!/bin/sh -e

#SBATCH --mem=4g
#SBATCH --output=Logs/3-reference/slurm-%A.out
#SBATCH --error=Logs/3-reference/slurm-%A.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}

# Document software versions used for publication
uname -a > Logs/3-reference/os-version-$SLURM_JOB_ID.txt 2>&1
samtools --version > Logs/3-reference/samtools-version-$SLURM_JOB_ID.txt 2>&1
gffread --version > Logs/3-reference/gffread-version-$SLURM_JOB_ID.txt 2>&1

#############################################################################
# Choose cdna.sh or gff2fasta.sh, or download the prebuilt kallisto
# file mus_musculus.tar.gz, extract, and copy transcriptome.idx to
# Data/3-kallisto-index/all-but-xy.index.
#
# There are multiple possible transcriptome references that can be used with
# kallisto.
#
# There are many gene IDs referenced in the GTF/GFF that are not in the CDNA.
# FIXME: Document the reason for this.  Are these predicted genes?
#
# There are a few gene IDs in the release 98 CDNA that are not in the GFF.
# If downstream analysis involved looking up genes in the GFF, this could
# result in a few misses.  This caused minor problems for CNC-EMDiff, which
# used CDNA as the reference.
#
# More importantly, the CDNA does not document features such as exons, UTRs,
# etc.  If downstream analysis will examine such features, use the GFF.

# Reference/cdna.sh       # Remove XY from Ensembl cdna transcriptome
Reference/gff2fasta.sh  # Construct genome minus XY using Ensembl GFF

Reference/create-chrom-sizes.sh
