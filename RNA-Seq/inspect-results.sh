#!/bin/sh -e

##########################################################################
#   Synopsis:
#       
#   Description:
#       Check logs and visualize results of each step in the pipeline
#
#   Arguments:
#       
#   Returns:
#
#   Examples:
#
#   Files:
#
#   Environment:
#
#   See also:
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
    if [ $# -lt 1 ]; then
	printf "Usage: view_files file-spec\n" >> /dev/stderr
	exit 1
    fi
    
    for file in $@; do
	printf "\n(V)iew $file\n(S)kip $file (default)\n(Q)uit\n"
	read view
	case $view in
	V|v)
	    more $file
	    ;;
	
	Q|q)
	    break
	    ;;
	esac
    done
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

wb_files()
{
    if [ $# -lt 1 ]; then
	printf "Usage: wb_files file-spec\n" >> /dev/stderr
	exit 1
    fi

    for file in $@; do
	printf "\n(V)iew $file?\n(S)kip $file (default)\n(Q)uit\n"
	read view
	case $view in
	V|v)
	    webbrowser $file
	    ;;
	
	Q|q)
	    break
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
Q.. Quit

EOM
    printf "Selection? "
    read selection
    
    case $selection in
    1)
	run_cmd "ls -l Data/01-organize/Raw-renamed | more"
	;;
    
    2)
	view_files Logs/02-qc-raw/*
	wb_files Data/02-qc-raw/*.html
	;;
    
    3)
	view_files Logs/03-multiqc-raw/*
	run_cmd "webbrowser ./Data/03-multiqc-raw/multiqc_report.html"
	;;

    4)
	view_files Logs/04-trim/*
	printf "Check for remaining adapter contamination? y/[n] "
	read scum
	if [ 0$scum = 0y ]; then
	    for fq in Data/04-trim/*.gz; do
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
	view_files Logs/05-qc-trimmed/*
	wb_files Data/05-qc-trimmed/*.html
	;;
    
    6)
	view_files Logs/06-multiqc-trimmed/*
	run_cmd "webbrowser ./Data/06-multiqc-trimmed/multiqc_report.html"
	;;

    7)
	view_files Logs/07-reference/*
	(cd Data/07-reference && view_files *.fa *.fai *.tsv)
	;;
    
    Q|q)
	;;
    
    *)
	printf "Invalid selection #: $selection\n" >> /dev/stderr
	;;

    esac
    pause
done

