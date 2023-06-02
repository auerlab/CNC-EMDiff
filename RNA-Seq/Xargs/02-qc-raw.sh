#!/bin/sh -e

# Raw files from sequencing center
raw=$1

# Filename stems for fastqc output
stem_raw=$(basename ${raw%.fastq.xz})
results=Results/02-qc-raw/$(basename ${stem_raw})_fastqc.html

if [ -e $results ]; then
    printf "$raw already processed.\n"
else
    printf "Processing $raw with fastqc...\n"
    xzcat $raw | fastqc -o Results/02-qc-raw stdin:$stem_raw \
	> Logs/02-qc-raw/$stem_raw.out 2> Logs/02-qc-raw/$stem_raw.err
fi
