#!/bin/sh -e

# Not yet sure how many threads are useful
# Maybe use --exclusive and all available cores?
#SBATCH --array=1-18 --nodes=1 --ntasks=2
# Andrea's script used 64000.  Just to force newer Mortimer nodes?
##SBATCH --mem=30000
# Based on 12 threads on peregrine
#SBATCH --mem=6g
##SBATCH --partition=batch,256g
#SBATCH --output=4-bwa-mem/slurm-%A_%a.out
#SBATCH --error=4-bwa-mem/slurm-%A_%a.err

# Fake SLURM_ARRAY_TASK_ID for testing script outside SLURM env
: ${SLURM_ARRAY_TASK_ID:=1}

cd 4-bwa-mem
infile1=../1-trim/*_S${SLURM_ARRAY_TASK_ID}_L001-R1.fastq.gz
infile2=$(echo $infile1 | sed -e 's|R1|R2|')
outfile=$(echo $infile1 | cut -d / -f 3 | sed -e 's|-R1.fastq.gz|.sam|')
ref_file=$(../Reference/ref-file-name.sh)
bwa mem -M -t ${SLURM_NTASKS} \
    ../Reference/$ref_file $infile1 $infile2 > $outfile
