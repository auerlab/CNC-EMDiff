#!/bin/sh -e

raw1=$1
base=$(basename $raw1)
stem=${base%%-R1.fastq.xz}
trimmed1=Results/04-trim/$stem-R1.fastq.gz

raw2=${raw1%%-R1.fastq.xz}-R2.fastq.xz
trimmed2=Results/04-trim/$stem-R2.fast1.gz

if [ -e $trimmed1 ]; then
    printf "$raw1 + $raw2 already processed.\n"
else
    printf "Processing $raw1 + $raw2 with fastq-trim...\n"
    
    log_stem=Logs/02-qc-raw/$stem_raw

    # Document software versions used for publication
    uname -a > $log_stem.out
    fastq-trim --version >> $log_stem.out
    pwd >> $log_stem.out

    if pwd | fgrep -q RNA-Seq; then
	adapter=AGATCGGAAGAG
    else
	adapter=CTGTCTCTTATACACATCT
    fi
    
    # fastq-trim is 2.5x faster with 1 core than cutadapt with 2 cores
    # Our reads use the default Illumina universal adapter,
    # but we'll state it explicitly anyway
    export GZIP=-1
    set -x
    time fastq-trim --3p-adapter1 $adapter --3p-adapter2 $adapter \
	--min-qual 24 --polya-min-length 4 $raw1 $trimmed1 $raw2 $trimmed2 \
	>> $log_stem.out 2>> $log_stem.err
fi
