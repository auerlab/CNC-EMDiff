#!/bin/sh -e

##########################################################################
#   Script description:
#       Trim adapters and low quality ends from Lumina reads
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#       All necessary tools are assumed to be in PATH.  If this is not
#       the case, add whatever code is needed here to gain access.
#       (Adding such code to your .bashrc or other startup script is
#       generally a bad idea since it's too complicated to support
#       every program with one environment.)
#
#   History:
#   Date        Name        Modification
#   2019-09-11  Jason Bacon Begin
##########################################################################

# Each job in the array will run a cutadapt (python) process and a
# compression process for part of the time.  If you don't want to
# oversubscribe compute nodes even for a little while, add --cpus-per-task=2
# There may be 2 pigz processes per job, but --cpus-per-task=3 doesn't help
# Limit concurrent jobs to avoid becoming I/O-bound.  Fastq-trim is so fast
# it ends up using only about 40% CPU while waiting for NFS on albacore.
# We'll finish the job array just as fast running only 9 at a time and
# getting 80% CPU utilization.
#SBATCH --array=1-18%9
# 2 xzcat, 2 gzip, and 1 fastq-trim, but xzcat and gzip only use about 50%
#SBATCH --cpus-per-task=2
# Memory requirements can only be determined by trial and error.
# Run a sample job and monitor closely in "top" or rununder a tool that
# reports maximum memory use.
#SBATCH --mem-per-cpu=250
#SBATCH --output=Logs/04-trim/slurm-%A_%a.out
#SBATCH --error=Logs/04-trim/slurm-%A_%a.err

# Set a default value for testing outside the SLURM environment
: ${SLURM_ARRAY_TASK_ID:=1}

# Document software versions used for publication
uname -a
# cutadapt --version
# fastq-trim --version
pwd

input1=$(ls Data/Raw-renamed/*-sample${SLURM_ARRAY_TASK_ID}-*R1*)
input2=$(ls Data/Raw-renamed/*-sample${SLURM_ARRAY_TASK_ID}-*R2*)
base=$(basename $input1)
stem=${base%%-R*.fastq.xz}

##############################################################################
# TrimGalore is just a wrapper around cutadapt and it cannot cut a fixed
# number of bases at the same time as trimming adapters ( -u and -U in
# cutadapt ), so we run cutadapt directly.
#
# trim_galore --illumina --stringency 3 -q 20 --paired -o Data/04-trim \
#    $input1 $input2 > Data/04-trim/${stem}-trash.out
# trim_galore --illumina --stringency 3 -q 20 --paired --trim1 -o Data/04-trim-trim1 \
#    $input1 $input2 > Data/04-trim-trim1/${stem}-trash.out
#
# Command from trim_galore log:
#   cutadapt -j 1 -e 0.1 -q 20 -O 3 -a AGATCGGAAGAGC
#   -j 1 -O 3 -e 0.1 are default for cutadapt

##############################################################################
# Remove ~15 bases at 5' end due to bias and last base at 3' end
# 5' bias is probably due to non-random cleavage and not a
# quality issue, so -u +15 -U +15 is probably not necessary
# -u +15 -U +15 \

# https://www.rootusers.com/gzip-vs-bzip2-vs-xz-performance-comparison/
# xz offers the best compression by far, but is slow at mid (-5) to high (-9)
# compression levels.  At -1, xz is faster than bzip2 while providing
# comparable compression.  If you want even faster compression and are
# willing to sacrifice compression ratio, use .gz or no compression for
# outputs instead.  A ZFS filesystem with lz4 compression should provide
# enough compression for intermediate files without gzip, bzip2, or xz.
# However, this may cause a network bottleneck as all processes write
# uncompressed FASTQ over NFS.  Using at least gzip -1 will transfer some
# of the load to the compute node CPUs and reduce NFS traffic considerably.
#
# suffix=.xz
suffix=.gz
output1=Data/04-trim/${stem}-R1.fastq$suffix
output2=Data/04-trim/${stem}-R2.fastq$suffix

# Maximize compression throughput so gzip is not a bottleneck.
# Don't use -1 if job is I/O bound.  Use the idle CPU to get better compression
# and reduce I/O.  The goal is to maximum CPU utilization of the fastq-trim
# process, as shown by job-top.  -4 seems to work best on albacore with its
# gigabit network.  Lower values are probably better with a high-speed network.
export GZIP=-4

if pwd | fgrep -q RNA-Seq; then
    adapter=AGATCGGAAGAG
else
    adapter=CTGTCTCTTATACACATCT
fi

# fastq-trim is 2.5x faster with 1 core than cutadapt with 2 cores
# Our reads use the default Illumina universal adatper,
# but we'll state it explicitly anyway
set -x
time fastq-trim --3p-adapter1 $adapter --3p-adapter2 $adapter \
    --min-qual 24 $input1 $output1 $input2 $output2

# To use cutadapt instead, comment out the above and uncomment this:
#cutadapt --compression-level=1
#    --quality-cutoff=24 \
#    -u -3 -U -3 \
#    -a $adapter -A $adapter \
#    --minimum-length=30 \
#    --output=Data/04-trim/${stem}-R1.fastq$suffix \
#    --paired-output=Data/04-trim/${stem}-R2.fastq$suffix \
#    $input1 $input2
