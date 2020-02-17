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

# Probably not much difference between these two methods
# Either works fine for a basic kallisto quantification
./cdna.sh           # Remove XY from Ensembl cdna transcriptome
# ./gtf2fasta.sh      # Construct genome-XY and use Ensembl GTF

reference=$(./reference-filename.sh)
samtools faidx $reference   # Needed for kiallisto --genomebam?

cd $save_cwd
kallisto index \
    --index=3-kallisto-index/all-but-xy.index \
    Reference/$reference
