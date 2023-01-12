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

view_logs()
{
    if [ $# != 1 ]; then
	printf "Usage: view_logs dir\n" >> /dev/stderr
	exit 1
    fi
    
    dir=$1
    for log in Logs/$dir/*; do
	more $log
	printf "View another? [y]/n "
	read another
	if [ 0$another = 0n ]; then
	    break;
	fi
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
4.. FastQ trim
Q.. Quit

EOM
    printf "Selection? "
    read selection
    
    case $selection in
    1)
	ls -l Data/01-organize/Raw-renamed | more
	;;
    
    2)
	view_logs 02-fastqc-raw
	for report in Data/02-qc-raw/*.html; do
	    webbrowser $report
	    printf "View another? [y]/n "
	    read another
	    if [ 0$another = 0n ]; then
		break;
	    fi
	done
	;;
    
    3)
	more Logs/03-multiqc-raw/*
	webbrowser ./Data/03-multiqc-raw/multiqc_report.html
	;;

    4)
	view_logs 04-trim
	printf "Check for remaining adapter contamination? y/[n] "
	read scum
	if [ 0$scum = 0y ]; then
	    for fq in Data/04-trim/*.gz; do
		$srun fastq-scum $fq
		printf "View another? [y]/n "
		read another
		if [ 0$another = 0n ]; then
		    break;
		fi
	    done
	fi
	;;
    
    *)
	printf "Invalid selection #: $selection\n" >> /dev/stderr
	exit 1
	;;
    esac
    pause
done

