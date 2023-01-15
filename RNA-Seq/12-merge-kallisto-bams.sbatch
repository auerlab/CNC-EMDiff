#!/bin/sh -e

#############################################################################
#   Description:
#       Merge BAM files for the same cell line and time point
#
#   Dependencies:
#       Requires kallisto pseudobams.  Run after *-kallisto-quant.sbatch.
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#############################################################################

# Set job array to number of samples / number of time points.
# 3 time points, chondro and neuro for each
#SBATCH --array=1-6
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
# top showed a steady 25M virtual, 12 res with a spike to 42/30 at the
# end (samtools index?) on FreeBSD 13
#SBATCH --mem=50m
#SBATCH --output=Logs/12-merge-kallisto-bams/slurm-%A_%a.out
#SBATCH --error=Logs/12-merge-kallisto-bams/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=4}

# Document software versions used for publication
uname -a
samtools --version
pwd

input_dir=Data/09-kallisto-quant
output_dir=Data/12-merge-kallisto-bams
if [ $SLURM_ARRAY_TASK_ID -le 3 ]; then
    cell_type=chondro
    time=$SLURM_ARRAY_TASK_ID
else
    cell_type=neuro
    time=$(($SLURM_ARRAY_TASK_ID - 3))
fi

input1=$(echo $input_dir/$cell_type-sample*-rep1-time$time/pseudoalignments.bam)
input2=$(echo $input_dir/$cell_type-sample*-rep2-time$time/pseudoalignments.bam)
input3=$(echo $input_dir/$cell_type-sample*-rep3-time$time/pseudoalignments.bam)

if [ $SLURM_ARRAY_TASK_ID -le 3 ]; then
    output=$output_dir/chondro-time$time.bam
else
    output=$output_dir/neuro-time$time.bam
fi

set -x
samtools merge -o $output $input1 $input2 $input3
samtools index $output
