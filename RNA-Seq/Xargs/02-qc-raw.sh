#!/bin/sh -e

# Raw files from sequencing center
raw=$1

printf "fastqc $raw...\n"

# Filename stems for fastqc output
stem_raw=$(basename ${raw%.fastq.xz})

xzcat $raw | fastqc -o Results/02-qc-raw stdin:$stem_raw

