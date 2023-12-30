#!/bin/sh -e

##########################################################################
#   Description:
#       Trim one raw fastq file.
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
raw1=$1
raw2=${raw1%%-R1.fastq.xz}-R2.fastq.xz

base=$(basename $raw1)
stem=${base%%-R1.fastq.xz}

# zstd is faster than gzip *and* provides a better compression ratio.
# Nobody should be using gzip anymore unless it's necessary for
# esoteric reasons.  Fastq-trim supports all common compression tools
# (gzip, bzip2, xz, zstd, lz4, ...).
trimmed1=Results/05-trim/$stem-R1.fastq.zst
trimmed2=Results/05-trim/$stem-R2.fastq.zst

if [ -e $trimmed1 ]; then
    printf "$raw1 + $raw2 already processed.\n"
else
    # Document software versions used for publication
    hostname
    uname -a
    fastq-trim --version
    pwd
    
    # fastq-trim is 2.5x faster with 1 core than cutadapt with 2 cores
    # Our reads use the default Illumina universal adapter,
    # but we'll state it explicitly anyway
    adapter=AGATCGGAAGAG
    
    set -x
    time fastq-trim --3p-adapter1 $adapter --3p-adapter2 $adapter \
	--min-qual 24 --polya-min-length 4 $raw1 $trimmed1 $raw2 $trimmed2
fi
