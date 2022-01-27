#!/bin/sh -e

##########################################################################
#   Remove identical reads from raw FASTQ files.  This will reduce the
#   amount of data processed by all subsequent steps and thus speed up
#   the analysis.
#
#   FIXME: Should this functionality be part of fastq-trim?
#
#   History:
#   Date        Name        Modification
#   2022-01-27  Jason Bacon Begin
##########################################################################

#SBATCH --array=1-18
#SBATCH --cpus-per-task=2
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem=2g
#SBATCH --output=Logs/02-fastq-derep/slurm-%A_%a.out
#SBATCH --error=Logs/02-fastq-derep/slurm-%A_%a.err

: ${SLURM_ARRAY_TASK_ID:=1}

infile1=$(ls Data/Raw-renamed/chondro-sample$SLURM_ARRAY_TASK_ID-*-R1.fastq.xz)
infile2=$(ls Data/Raw-renamed/chondro-sample$SLURM_ARRAY_TASK_ID-*-R2.fastq.xz)

outfile1=Data/02-fastq-derep/$(basename ${infile1%.xz}).gz
outfile2=Data/02-fastq-derep/$(basename ${infile2%.xz}).gz

printf "Dereplicating $infile1 -> $outfile1...\n"
xzcat $infile1 | blt fastx-derep | gzip -1 > $outfile1

printf "Dereplicating $infile1 -> $outfile1...\n"
xzcat $infile2 | blt fastx-derep | gzip -1 > $outfile2
