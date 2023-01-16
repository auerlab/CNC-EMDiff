#!/bin/sh -e

##########################################################################
#   Description:
#       Help the user identify biologically significant fold-changes
#       in the DE output.  This involves filtering based on
#       user-provided criteria, showing fold-change, P-value, other
#       statistics, available gene function information, etc.
#
#   History:
#   Date        Name        Modification
#   2023-01-12  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 \n"
    exit 1
}


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


##########################################################################
#   Main
##########################################################################

if [ $# != 0 ]; then
    usage
fi

# Possibly use GAF files directly, or tools sets like goatools or ermineJ

gaf_file=goa_human.gaf
obo_file=go.obo

if [ ! -e $gaf_file ]; then
    curl -O http://current.geneontology.org/annotations/$gaf_file.gz
    gunzip $gaf_file.gz
fi

if [ ! -e $obo_file ]; then
    curl -O http://current.geneontology.org/ontology/go.obo
fi

awk -f go-names.awk go.obo > go-names.tsv

while true; do
    clear
    cat << EOM

1.. Browse Kallisto DE results
2.. Filter Kallisto DE results
Q.. Quit

EOM
    printf "Selection? "
    read selection
    
    case 0$selection in
    01)
	more Data/13-fasda-kallisto/*.txt
	;;
    
    02)
	cat << EOM
			    **** NOTE ****
			    
What constitutes a meaningful read count varies immensely from one gene
to another.  It also depends directly on the count of DNA/RNA fragments
in the sample, which does not reflect the count per cell.  Minimum high
read count is a arbitrary decision unless you know the parameters of
your read data and genes of interest.  Use with caution.

What is a meaningful fold-change also varies immensely among genes.
1.5 might be meaningful for an enzyme whose affect is amplified, and
meaningless for genes whose translation level is limited by other factors.

P-values represent STATISTICAL significance of a fold-change based on
unreliable read counts, not BIOLOGICAL significance of a change in the cell.
Read counts are highly unreliable due to the many uncontrollable variables
in the cell and the sequencing process that distort read counts.

Users should not rely on any one set of parameters here to determine
which fold-changes are biologically significant.

EOM
	printf "Minimum high read count? [30] "
	read min_count
	if [ -z "$min_count" ]; then
	    min_count=30
	fi

	printf "Minimum %% agreement across replicates for up/down regulation? [66] "
	read min_agree
	if [ -z "$min_agree" ]; then
	    min_agree=66
	fi
	
	printf "Mininum FC? [2] "
	read min_fc
	if [ -z "$min_fc" ]; then
	    min_fc=2
	fi
	
	printf "Maximum P-value? [0.2] "
	read max_pv
	if [ -z "$max_pv" ]; then
	    max_pv=0.2
	fi
	
	# Need to map gene IDs to UniProt IDs in order to search GO
	# associations
	# printf "GO term? "
	# read go_term

	printf "%-30s %-s\n" "File" "Features meeting criteria"
	for file in Data/13-fasda-kallisto/*.txt; do
	    printf "%-30s" `basename $file`
	    awk -v min_count=$min_count \
		-v min_agree=$min_agree \
		-v min_fc=$min_fc \
		-v max_pv=$max_pv \
		'(($2 >= min_count) || ($3 >= min_count)) && ($6 >= min_agree) && ($7 >= min_fc) && ($8 <= max_pv)' \
		$file | wc -l
	done
	printf "\nView results? [y]/n "
	read view
	if [ 0$view != 0n ]; then
	    for file in Data/13-fasda-kallisto/*.txt; do
		basename $file
		awk -v min_count=$min_count \
		    -v min_agree=$min_agree \
		    -v min_fc=$min_fc \
		    -v max_pv=$max_pv \
		    '(($2 >= min_count) || ($3 >= min_count)) && ($6 >= min_agree) && ($7 >= min_fc) && ($8 <= max_pv)' \
		    $file | more
	    done
	fi
	;;
    
    0Q|0q)
	exit 0
	;;
    
    *)
	printf "Invalid selection.\n"
	pause
    
    esac
done

