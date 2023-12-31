#!/bin/sh -e

##########################################################################
#   Script description:
#       Run quality checks on one raw data file.
##########################################################################

usage()
{
    printf "Usage: $0 raw-file.fastq\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# -ne 1 ]; then
    usage
fi
fastq=$1

# Filename stems for fastqc output
stem_fastq=$(basename ${fastq%.fastq.xz})
results=Results/03-qc-raw/$(basename ${stem_fastq})_fastqc.html

if [ -e $results ]; then
    printf "$fastq already processed.\n"
else
    printf "Processing $fastq with fastqc...\n"
    
    # Document software versions used for publication
    hostname
    uname -a
    fastqc --version
    pwd
    date
    xzcat $fastq | fastqc -o Results/03-qc-raw stdin:$stem_fastq
    date
fi
