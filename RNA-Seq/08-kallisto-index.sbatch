#!/bin/sh -e

##########################################################################
#   Description:
#       Build kallisto index for reference transcriptome.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#       
#   History:
#   Date        Name        Modification
#   2019-11-??  Jason Bacon Begin
##########################################################################

# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem=4g
#SBATCH --output=Logs/08-kallisto-index/slurm-%A.out
#SBATCH --error=Logs/08-kallisto-index/slurm-%A.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
kallisto version
samtools --version
pwd

transcriptome=Data/07-reference/$(Reference/transcriptome-filename.sh)
printf "Using reference $transcriptome...\n"

# Needed for kallisto --genomebam
if [ ! -e $transcriptome.fai ]; then
    printf "Building $transcriptome...\n"
    samtools faidx $transcriptome
fi

printf "Building kallisto index...\n"
set -x
kallisto index --index=Data/08-kallisto-index/all-but-xy.index $transcriptome
