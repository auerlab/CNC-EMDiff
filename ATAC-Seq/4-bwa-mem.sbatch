#!/bin/sh -e

# Not yet sure how many threads are useful
# Maybe use --exclusive and all available cores?
# If cores are limited, smaller thread counts are usually more efficient
# as they will require less communication overhead
#SBATCH --array=1-18 --cpus-per-task=2
# Andrea's script used 64000.  Just to force newer Mortimer nodes?
##SBATCH --mem=30000
# Based on 2 threads on albacore
# 12 threads on peregrine used close to 6g
#SBATCH --mem=5g
##SBATCH --partition=batch,256g 
# --exclude=compute-[001-008,012]
#SBATCH --output=4-bwa-mem/slurm-%A_%a.out
#SBATCH --error=4-bwa-mem/slurm-%A_%a.err

# Spoof SLURM_ARRAY_TASK_ID if not set for testing outside SLURM env
: ${SLURM_ARRAY_TASK_ID:=1}

cd 4-bwa-mem
infile1=../1-trim/*_S${SLURM_ARRAY_TASK_ID}_L001-R1.fastq.gz
infile2=../1-trim/*_S${SLURM_ARRAY_TASK_ID}_L001-R2.fastq.gz
outfile=$(echo $infile1 | cut -d / -f 3 | sed -e 's|-R1.fastq.gz|.sam|')
ref_file=$(../Reference/ref-file-name.sh)

srun bwa mem -M -t $SLURM_CPUS_PER_TASK \
    ../Reference/$ref_file $infile1 $infile2 > $outfile
