#!/bin/sh -e

# Must be run from parent parent of Verify dir for relative paths to work
if [ `basename $(pwd)` == Verify ]; then
    cd ..
fi

ref_dir=Results/07-reference
gtf=$(Reference/gtf-filename.sh)

printf "Transcripts in GTF file:\n"
awk '$3 ~ "transcript"' $ref_dir/$gtf | wc -l

printf "\nThere should be 6 normalized abundance files and 6 FC files.\n"

printf "\nLine counts in normalized abundance files:\n"
printf "Should be the same as above.\n\n"
wc -l Results/13-fasda-kallisto/*-all-norm.tsv

printf "\nLine counts in fold-change files:\n"
printf "Should be one more than above (to account for header line).\n\n"
wc -l Results/13-fasda-kallisto/*-FC.txt

# Add more sophisticated checks here if desired
