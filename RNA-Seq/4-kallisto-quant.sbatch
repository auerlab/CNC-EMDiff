#!/bin/sh -e

#SBATCH --mem-per-cpu=4g
#SBATCH --array=1-18
#SBATCH --output=4-kallisto-quant/slurm-%A_%a.out
#SBATCH --error=4-kallisto-quant/slurm-%A_%a.err
#SBATCH --exclude=compute-[001-008,012]
##SBATCH --partition=256g

# CentOS pkgsrc
if [ -e /sharedapps/pkg-2019Q2/bin/kallisto ]; then
    # sh is bash on CentOS, so this is OK
    source /etc/bashrc
    module load /sharedapps/pkg-2019Q2/etc/modulefiles/pkgsrc/2019Q2
fi

# Run kallisto with 500 bootstraps for Sleuth
mkdir -p 4-kallisto-quant/$SLURM_ARRAY_TASK_ID
kallisto quant \
    -i 3-kallisto-index/all-but-xy.index \
    -o 4-kallisto-quant/$SLURM_ARRAY_TASK_ID \
    -b 500 \
    1-trim/*_S${SLURM_ARRAY_TASK_ID}_L002-R1.fastq.gz \
    1-trim/*_S${SLURM_ARRAY_TASK_ID}_L002-R2.fastq.gz
