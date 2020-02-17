#!/bin/sh -e

# Uses all available RAM
##SBATCH --exclusive
#SBATCH --mem=4g
##SBATCH --partition=256g
#SBATCH --output=3-kallisto-index/slurm-%A.out
#SBATCH --error=3-kallisto-index/slurm-%A.err

# CentOS pkgsrc
if [ -e /sharedapps/pkg-2019Q2/bin/kallisto ]; then
    # sh is bash on CentOS, so this is OK
    source /etc/bashrc
    module load /sharedapps/pkg-2019Q2/etc/modulefiles/pkgsrc/2019Q2
fi

save_cwd=$(pwd)
cd Reference

############################################################################
# Choose cdna or gtf2fasta or copy the kallisto transcriptome.idx file from
# ../Reference/mus_musculus.tar.gz to 3-kallisto-index/all-but-xy.index
# and skip running this sbatch script.

./cdna.sh           # Remove XY from Ensembl cdna transcriptome
# ./gtf2fasta.sh      # Construct genome-XY and use Ensembl GTF

reference=$(./reference-filename.sh)
samtools faidx $reference   # Needed for kiallisto --genombam?

cd $save_cwd
kallisto index \
    --index=3-kallisto-index/all-but-xy.index \
    Reference/$reference
