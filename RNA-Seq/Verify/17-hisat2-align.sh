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

more Logs/17-hisat2-align/*.out
more Logs/17-hisat2-align/*.err

# Add more sophisticated checks here if desired
