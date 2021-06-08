#!/bin/sh -e

#SBATCH --mem=4g
#SBATCH --output=Data/3-kallisto-index/slurm-%A.out
#SBATCH --error=Data/3-kallisto-index/slurm-%A.err

save_cwd=$(pwd)
cd Reference

#############################################################################
# Choose cdna or gtf2fasta, or copy the kallisto transcriptome.idx file from
# ../Reference/mus_musculus.tar.gz to Data/3-kallisto-index/all-but-xy.index
# and skip running this sbatch script.

# There are multiple possible transcriptome references that can be used with
# kallisto.
#
# There are many gene IDs references in the GTF/GFF that are not in the CDNA.
# FIXME: Document the reason for this.  Are these predicted genes?
#
# There are a few gene IDs in the release 98 CDNA that are not in the GTF.
# If downstream analysis involved looking up genes in the GTF, this could
# result in a few misses.  This caused minor problems for CNC-EMDiff, which
# used CDNA as the reference.
#
# More importantly, the CDNA does not document features such as exons, UTRs,
# etc.  If downstream analysis will examine such features, use the GTF.
#

./cdna.sh           # Remove XY from Ensembl cdna transcriptome
# ./gtf2fasta.sh      # Construct genome minus XY using Ensembl GTF

reference=$(./reference-filename.sh)
samtools faidx $reference   # Needed for kiallisto --genomebam?

cd $save_cwd
kallisto index \
    --index=Data/3-kallisto-index/all-but-xy.index \
    Reference/$reference
