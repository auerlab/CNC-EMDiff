#!/bin/sh -e

#SBATCH --mem=4g
#SBATCH --output=3-kallisto-index/slurm-%A.out
#SBATCH --error=3-kallisto-index/slurm-%A.err

# Mortimer pkgsrc
if [ -e /sharedapps/pkg-2019Q2/bin/kallisto ]; then
    # sh is bash on CentOS, so this is OK
    source /etc/bashrc
    module load /sharedapps/pkg-2019Q2/etc/modulefiles/pkgsrc/2019Q2
fi

save_cwd=$(pwd)
cd Reference

#############################################################################
# Choose cdna or gtf2fasta, or copy the kallisto transcriptome.idx file from
# ../Reference/mus_musculus.tar.gz to 3-kallisto-index/all-but-xy.index
# and skip running this sbatch script.

# Not much difference between these two methods. Either works fine for a basic
# kallisto quantificatiod. There are a few gene IDs in the release 98 CDNA
# reference that are not in the GTF.  If downstream analysis involved looking
# up genes in the GTF, this could result in a few misses.
./cdna.sh           # Remove XY from Ensembl cdna transcriptome
# ./gtf2fasta.sh      # Construct genome minus XY using Ensembl GTF

reference=$(./reference-filename.sh)
samtools faidx $reference   # Needed for kiallisto --genomebam?

cd $save_cwd
kallisto index \
    --index=3-kallisto-index/all-but-xy.index \
    Reference/$reference
