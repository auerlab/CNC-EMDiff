#!/bin/sh -e

##########################################################################
#   Synopsis:
#       
#   Description:
#       Check logs and visualize results of each step in the pipeline
#
#   History:
#   Date        Name        Modification
#   2023-01-12  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0\n"
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
#   Function description:
#       
#   Arguments:
#       
#   Returns:
#       
#   History:
#   Date        Name        Modification
#   2023-01-12  Jason Bacon Begin
##########################################################################

run_cmd()
{
    if [ $# != 1 ]; then
	printf 'Usage: run_cmd "command"\n' >> /dev/stderr
	exit 1
    fi
    
    printf "Running '%s':\n" "$1"
    pause
    eval $1
    return 0
}


##########################################################################
#   Function description:
#       
#   Arguments:
#       
#   Returns:
#       
#   History:
#   Date        Name        Modification
#   2023-01-12  Jason Bacon Begin
##########################################################################

view_files()
{
    if [ $# -lt 2 ]; then
	printf 'Usage: view_files "command" file-spec\n' >> /dev/stderr
	exit 1
    fi
    command="$1"
    shift
    
    for file in $@; do
	printf "\n(V)iew $file (default)\n(S)kip\n(Q)uit\n"
	read view
	case $view in
	S|s)
	    ;;
	
	Q|q)
	    break
	    ;;
	
	*)
	    $command $file
	    ;;
	esac
    done
    return 0
}


##########################################################################
#   Main
##########################################################################

if [ $# != 0 ]; then
    usage
fi

if which srun; then
    srun=srun
fi

selection=''
while [ 0$selection != 0q ]; do
    clear
    cat << EOM

			=================================
			    RNA-Seq results inspector
			=================================
			
1.. Raw renamed links
2.. Raw FastQC reports
3.. Raw MultiQC report
4.. Adapter trimming
5.. Trimmed FastQC reports
6.. Trimmed MultiQC report
7.. Reference transcriptome
8.. Kallisto index
9.. Kallisto abundances
10. SAM FastQC reports
11. SAM MultiQC report
12. Merged kallisto BAMs (IGV)
13. FASDA DE results for kallisto
16. Hisat2 index
17. Hisat2 alignments
19. FASDA abundance for hisat2
20. FASDA DE results for hisat2 
Q.. Quit

EOM
    printf "Selection? "
    read selection
    
    case $selection in
    1)
	run_cmd "ls -l Results/01-organize/Raw-renamed | more"
	;;
    
    2)
	view_files more Logs/02-qc-raw/*
	view_files webbrowser Results/02-qc-raw/*.html
	;;
    
    3)
	view_files more Logs/03-multiqc-raw/*
	run_cmd "webbrowser ./Results/03-multiqc-raw/multiqc_report.html"
	;;

    4)
	view_files more Logs/04-trim/*
	printf "Check for remaining adapter contamination (may take a long time)? y/[n] "
	read scum
	if [ 0$scum = 0y ]; then
	    for fq in Results/04-trim/*.gz; do
		run_cmd "$srun fastq-scum $fq"
		printf "View another? [y]/n "
		read another
		if [ 0$another = 0n ]; then
		    break;
		fi
	    done
	fi
	;;
    
    5)
	view_files more Logs/05-qc-trimmed/*
	view_files webbrowser Results/05-qc-trimmed/*.html
	;;
    
    6)
	view_files more Logs/06-multiqc-trimmed/*
	run_cmd "webbrowser ./Results/06-multiqc-trimmed/multiqc_report.html"
	;;

    7)
	view_files more Logs/07-reference/*
	(cd Results/07-reference && view_files more *.fa *.fai *.tsv)
	;;
    
    8)
	view_files more Logs/08-kallisto-index/*
	run_cmd "$srun kallisto inspect Results/08-kallisto-index/all-but-xy.index | more"
	;;
    
    9)
	view_files more Logs/09-kallisto-quant/*
	view_files more Results/09-kallisto-quant/*/abundance.tsv
	;;

    10)
	view_files more Logs/10-qc-sam/*
	view_files webbrowser Results/10-qc-sa/*.html
	;;
    
    11)
	view_files more Logs/11-multiqc-sam/*
	run_cmd "webbrowser ./Results/11-multiqc-sam/multiqc_report.html"
	;;

    12)
	view_files more Logs/12-merge-kallisto-bams/*
	if which igv > /dev/null 2>&1; then
	    (cd Results/12-merge-kallisto-bams && igv)
	else
	    printf "IGV is not installed on $(hostname).\n"
	    printf "Please install IGV and try again.\n"
	fi
	;;
    
    13)
	view_files more Logs/13-fasda-kallisto/*
	view_files more Results/13-fasda-kallisto/*.txt
	;;

    16)
	view_files more Logs/16-hisat2-index/*
	# more doesn't function under srun, so use cat and pipe
	# through local more process
	run_cmd "$srun cat Results/16-hisat2-index/all-but-xy.genome.fa | more"
	run_cmd "$srun cat Results/16-hisat2-index/all-but-xy.genome.fa.fai | more"
	;;
    
    17)
	view_files more Logs/17-hisat2-align/*
	# FIXME: srun has trouble with "more"
	#for file in Results/17-hisat2-align/*.bam; do
	#    run_cmd "$srun samtools view $file" | more
	#done
	;;
	
    19)
	view_files more Logs/19-fasda-abundance-hisat2/*
	# FIXME: Move update FASDA to place abundance files in separate dir
	view_files more Results/17-hisat2-align/*.tsv
	;;

    20)
	view_files more Logs/20-fasda-fc-hisat2/*
	view_files more Results/20-fasda-fc-hisat2/*.txt
	;;

    Q|q)
	exit 0
	;;
    
    *)
	printf "Invalid selection #: $selection\n" >> /dev/stderr
	;;

    esac
    pause
done

