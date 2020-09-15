#!/bin/sh -e

ref_dir=../../RNA-Seq/Reference
gtf=$ref_dir/Mus_musculus.GRCm38.98.gtf
awk '{ print $3 }' $gtf | sort | uniq
