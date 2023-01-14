#!/bin/sh -e

##########################################################################
#   Description:
#       Run STAR aligner on each RNA sample.
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

# Set job array to number of samples.
#SBATCH --array=1-18
# 1 for STAR, 1 for each xzcat
#SBATCH --cpus-per-task=3
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem=24000
#SBATCH --output=Logs/21-star-align/slurm-%A_%a.out
#SBATCH --error=Logs/21-star-align/slurm-%A_%a.err

##############################################################################
# Align with star, which can handle splice junctions in RNA reads

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_CPUS_PER_TASK:=4}

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

build=$(../Common/genome-build.sh)
release=$(../Common/genome-release.sh)
genome=$(Reference/genome-filename.sh)

xz1=$(echo Data/03-trim/*sample${SLURM_ARRAY_TASK_ID}-*-R1.fastq.xz)
xz2=$(echo Data/03-trim/*sample${SLURM_ARRAY_TASK_ID}-*-R2.fastq.xz)

xzb=$(basename $xz1)
filename_prefix=${xzb%-R*}-
bam=Data/21-star-align/${filename_prefix}Aligned.sortedByCoord.out.bam

# Empty BAM files may be caused by failing --readFilesCommand
# This didn't help.  2 issues on Github about empty BAM files, so
# probably a bug.
# test -e $unxz1 || xzcat $xz1 > $unxz1
# test -e $unxz2 || xzcat $xz2 > $unxz2

STAR \
    --genomeDir Data/20-star-index \
    --runThreadN 1 \
    --readFilesIn $xz1 $xz2 \
    --readFilesCommand xzcat \
    --outFileNamePrefix Data/21-star-align/$filename_prefix \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --outSAMattributes Standard

#samtools index $bam
