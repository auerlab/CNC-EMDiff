#!/bin/sh -e

# Not yet sure how many threads are useful
#SBATCH --array=1-18 --nodes=1 --ntasks=12
# Andrea's script used 64000.  Just to force newer Mortimer nodes?
##SBATCH --mem=30000
# Based on 12 threads on peregrine
#SBATCH --mem=6g
#SBATCH --partition=batch,256g

# Fake SLURM_ARRAY_TASK_ID for testing script outside SLURM env
: ${SLURM_ARRAY_TASK_ID:=1}

infile1=1-trim/*_S${SLURM_ARRAY_TASK_ID}_L001-R1.fastq.gz
infile2=$(echo $infile1 | sed -e 's|R1|R2|')
outfile=4-BWA-mem/$(echo $infile1 | cut -d / -f 2 | sed -e 's|-R1.fastq.gz|.sam|')
ref_file=$(./ref-file-name)

bwa mem -M -t ${SLURM_NTASKS} $ref_file $infile1 $infile2 > $outfile
