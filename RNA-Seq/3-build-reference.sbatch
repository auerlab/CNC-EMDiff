#!/bin/sh -e

#SBATCH --mem=4g
#SBATCH --output=Logs/3-reference/slurm-%A.out
#SBATCH --error=Logs/3-reference/slurm-%A.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}

# Document software versions used for publication
uname -a > Logs/3-reference/os-version-$SLURM_JOB_ID.txt 2>&1

#############################################################################
# Choose cdna or gtf2fasta, or copy the kallisto transcriptome.idx file from
# ../Reference/mus_musculus.tar.gz to Data/3-kallisto-index/all-but-xy.index
# and skip running this sbatch script.
#
# There are multiple possible transcriptome references that can be used with
# kallisto.
#
# There are many gene IDs referenced in the GTF/GFF that are not in the CDNA.
# FIXME: Document the reason for this.  Are these predicted genes?
#
# There are a few gene IDs in the release 98 CDNA that are not in the GTF.
# If downstream analysis involved looking up genes in the GTF, this could
# result in a few misses.  This caused minor problems for CNC-EMDiff, which
# used CDNA as the reference.
#
# More importantly, the CDNA does not document features such as exons, UTRs,
# etc.  If downstream analysis will examine such features, use the GTF.

# Reference/cdna.sh           # Remove XY from Ensembl cdna transcriptome
Reference/gtf2fasta.sh      # Construct genome minus XY using Ensembl GTF

Reference/create-chrom-sizes.sh
