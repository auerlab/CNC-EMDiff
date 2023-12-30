#!/bin/sh -e

##########################################################################
#   Example verification script.  Edit to add additional checks or
#   tune the existing checks to your personal taste.
##########################################################################

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}

cmd="webbrowser Results/04-multiqc-raw/multiqc_report.html"
printf "Running '$cmd'...\n"
pause
$cmd
