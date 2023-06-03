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
    
    log_stem=Logs/02-qc-raw/$stem_raw
    # Document software versions used for publication
    uname -a > $log_stem.out
    fastqc --version >> $log_stem.out
    pwd $log_stem.out

    xzcat $raw | fastqc -o Results/02-qc-raw stdin:$stem_raw \
	>> $log_stem.out 2>> $log_stem.err
fi
