#!/bin/sh -e


##########################################################################
#   Function description:
#       Pause until user presses return
##########################################################################

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}

#export PATH=$HOME/Prog/Src/local/bin:$PATH
which fasda

dir=Data/09-kallisto-quant
for cell_type in neuro chondro; do
    for condition in time1 time2; do
	printf "Normalizing $condition...\n"
	time fasda normalize --output $cell_type-$condition-all-norm.tsv \
	    $dir/$cell_type-*-rep*-$condition/abundance.tsv
    done
    
    printf "Computing fold-change...\n"
    time fasda fold-change --output $cell_type-time1-time2-FC.txt \
	$cell_type-time1-all-norm.tsv $cell_type-time2-all-norm.tsv
    pause
    more $cell_type-time1-time2-FC.txt
done
