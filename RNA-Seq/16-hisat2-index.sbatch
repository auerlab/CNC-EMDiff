#!/bin/sh -e

##########################################################################
#   Description:
#       Build hisat2 index for reference genome.
#
#   Dependencies:
#       Requires a reference genome.  Run after *-reference.sbatch.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#       
#   History:
#   Date        Name        Modification
#   2021-11-24  Jason Bacon Begin
##########################################################################

#SBATCH --mem=6g
#SBATCH --cpus-per-task=4
#SBATCH --output=Logs/16-hisat2-index/slurm-%A.out
#SBATCH --error=Logs/16-hisat2-index/slurm-%A.err

# Document software versions used for publication
uname -a
hisat2 --version
samtools --version
pwd

# Run hisat2-build on a copy in 16-hisat2-index so it will put the .ht2
# files there
genome=$(Reference/genome-filename.sh)
ln -f Results/07-reference/$genome Results/16-hisat2-index
genome=Results/16-hisat2-index/$genome
printf "Using reference $genome...\n"

# FIXME: Do this in *-reference.sbatch
if [ ! -e $genome.8.ht2 ]; then
    printf "Building $genome.*.ht2...\n"
    set -x
    hisat2-build -p $SLURM_CPUS_PER_TASK $genome $genome
    set +x
fi
if [ ! -e $genome.fai ]; then
    printf "Building $genome.fai...\n"
    set -x
    samtools faidx $genome
    set +x
fi
ls Results/16-hisat2-index
