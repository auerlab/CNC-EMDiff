#!/bin/sh -e

# Must be run from parent parent of Verify dir for relative paths to work
if [ `basename $(pwd)` == Verify ]; then
    cd ..
fi

ref_dir=Results/07-reference
gtf=$(Reference/gtf-filename.sh)

printf "Transcripts in GTF file:\n"
awk '$3 ~ "transcript"' $ref_dir/$gtf | wc -l

printf "\nLine counts in kallisto abundance files:\n"
printf "Should be 1 more than above (to account for header line).\n\n"
wc -l Results/09-kallisto-quant/*/abundance.tsv
