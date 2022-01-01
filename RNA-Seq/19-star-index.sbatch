#!/bin/sh -e

##########################################################################
#   Description:
#       Build STAR aligner index for reference genome.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#       
#   History:
#   Date        Name        Modification
#   2021-11-26  Jason Bacon Begin
##########################################################################

##########################################################################
#   https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html
##########################################################################

#SBATCH --cpus-per-task=6
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem=60g
#SBATCH --output=Logs/19-star-index/slurm-%A.out
#SBATCH --error=Logs/19-star-index/slurm-%A.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_CPUS_PER_TASK:=6}

cat << EOM

STAR genomeGenerate and alignment both report success, but all BAM output
files are empty.  This is apparently a common problem that the developers
are aware of but do not know how to solve.

In addition, STAR requires more than 30 gigabytes of memory for some
operations.

Use of hisat2 is recommended at this time.  It is faster, more reliable,
and more memory efficient.

EOM
exit

# Document software versions used for publication
uname -a
STAR --version
pwd

genome=$(Reference/genome-filename.sh)
gtf=$(Reference/gtf-filename.sh)
printf "Using reference $genome, $gtf...\n"

# Overhang should be read length - 1
STAR --runThreadN $SLURM_CPUS_PER_TASK \
    --runMode genomeGenerate \
    --genomeDir Data/19-star-index \
    --genomeFastaFiles Data/03-reference/$genome \
    --sjdbGTFfile Data/03-reference/$gtf \
    --sjdbOverhang 97

ls Data/19-star-index
