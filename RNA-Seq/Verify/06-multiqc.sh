#!/bin/sh -e

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}

# Must be run from parent parent of Verify dir for relative paths to work
if [ `basename $(pwd)` == Verify ]; then
    cd ..
fi

firefox Results/06-multiqc-trimmed/multiqc_report.html

# Add more sophisticated checks here if desired
