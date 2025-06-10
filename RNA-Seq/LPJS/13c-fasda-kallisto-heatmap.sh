#!/bin/sh -e

##########################################################################
#   Description:
#       Generate and display heatmap for Kallisto results
#       
#   History:
#   Date        Name        Modification
#   2025-03-31  Jason Bacon Begin
##########################################################################

input_dir=Results/09-kallisto-quant
samples=$(ls -d $input_dir/chondro* | wc -l)
replicates=$(($samples / 3))

for cell_type in chondro neuro; do
    for start_time in $(seq 1 2); do
	for end_time in $(seq $((start_time + 1)) 3); do
	    printf "Comparing time$start_time and time$end_time...\n"
	    # Generate a list of features of interest by filtering the FASDA fold-changes
	    features=filtered-features.txt
	    # chondro time2-time3 has nothing with P < 0.05 (odd)
	    fasda filter --max-p-val 0.10 \
		Results/13-fasda-kallisto/$cell_type-time$start_time-time$end_time-FC.txt \
		| awk '{ print $1 }' | head -n 30 > $features
	    wc -l $features
	    
	    # FIXME: Update when ./heatmap is moved into PATH
	    # Add --debug to see Python data structures
	    fasda heatmap $features \
		Results/13-fasda-kallisto/$cell_type-time$start_time-norm.tsv \
		Results/13-fasda-kallisto/$cell_type-time$end_time-norm.tsv
	done
    done
done
