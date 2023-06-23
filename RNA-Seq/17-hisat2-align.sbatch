#!/bin/sh -e

##########################################################################
#   Description:
#       Run hisat2 aligner on each RNA sample.
#
#   Dependencies:
#       Requires hisat2 index.  Run after *-hisat2-index.sbatch.
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

# Set job array to number of samples.
#SBATCH --array=1-18
#SBATCH --cpus-per-task=4
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or run under a tool that
# reports maximum memory use.
#SBATCH --mem=4000
#SBATCH --output=Logs/17-hisat2-align/slurm-%A_%a.out
#SBATCH --error=Logs/17-hisat2-align/slurm-%A_%a.err

##############################################################################
# Align with hisat2, which can handle splice junctions in RNA reads

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
hisat2 --version
pwd

genome=$(Reference/genome-filename.sh)

# samtools sort dumps temp files in CWD
cd Results/17-hisat2-align

gz1=$(echo ../04-trim/*sample${SLURM_ARRAY_TASK_ID}-*-R1.fastq.gz)
gz2=$(echo ../04-trim/*sample${SLURM_ARRAY_TASK_ID}-*-R2.fastq.gz)

gzb=$(basename $gz1)
bam=${gzb%-R*}.bam

set -x
hisat2 --threads $SLURM_CPUS_PER_TASK -x ../16-hisat2-index/$genome \
    -1 $gz1 -2 $gz2 | samtools sort > $bam
samtools index $bam
