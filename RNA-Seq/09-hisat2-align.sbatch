#!/bin/sh -e

# FIXME: mem is from kallisto
#SBATCH --mem=6000
#SBATCH --array=1-18
#SBATCH --ntasks=2
#SBATCH --output=Logs/09-hisat2-align/slurm-%A_%a.out
#SBATCH --error=Logs/09-hisat2-align/slurm-%A_%a.err

##############################################################################
# Align with hisat2, which can handle splice junctions in RNA reads

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}
: ${SLURM_JOB_ID:=test}

# Document software versions used for publication
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
    uname -a > Logs/09-hisat2-align/os-version-$SLURM_JOB_ID.txt 2>&1
    hisat2 --version > Logs/09-hisat2-align/hisat2-version-$SLURM_JOB_ID.txt 2>&1
fi

build=$(../Common/genome-build.sh)
release=$(../Common/genome-release.sh)
genome=$(Reference/genome-filename.sh)

cd Data/09-hisat2-align

# hisat2 can't handle xz.
# Convert xz to gz rather than raw to reduce NFS load from compute nodes
xz1=$(echo ../01-trim/*-sample${SLURM_ARRAY_TASK_ID}-*-R1.fastq.xz)
xz2=$(echo ../01-trim/*-sample${SLURM_ARRAY_TASK_ID}-*-R2.fastq.xz)
gz1=${xz1%.xz}.gz
gz2=${xz2%.xz}.gz
test -e $gz1 || xzcat $xz1 | gzip --stdout > $gz1
test -e $gz2 || xzcat $xz2 | gzip --stdout > $gz2

bam=$(basename $gz1)
bam=${bam%-R*}.bam

set -x
hisat2 --threads 2 -x ../../Data/hisat2-index/$genome \
    -1 $gz1 -2 $gz2 | samtools sort > $bam
samtools index $bam

# Remove temporary .gz files, keep .xz which is ~40% smaller
rm -f $gz1 $gz2
