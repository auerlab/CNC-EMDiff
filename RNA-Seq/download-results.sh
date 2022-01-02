#!/bin/sh -e

##########################################################################
#   Description:
#       Display rsync commands to copy and paste for downloading results
#       
#   History:
#   Date        Name        Modification
#   2022-01-02  Jason Bacon Begin
##########################################################################

usage()
{
    cat << EOM

Usage: $0 remote-hostname

remote-hostname is either the external hostname of this machine (for
pulling) or the hostname of the destination (for pushing).

EOM
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

if [ $# != 1 ]; then
    usage
fi
remote_host=$1

if [ $0 != ./download-results.sh ]; then
    printf "$0 should be run as ./download-results.sh\n"
    exit 1
fi

cat << EOM

Copy and paste the commands below to pull or push results from $(hostname -s).

EOM
pause
scripts=$(ls 0[2-9]-* [1-9][0-9]-*)
for script in $scripts; do
    stage=${script%.*}
    cat << EOM
---
rsync -av $remote_host:$(pwd)/Data/$stage .
rsync -av $(pwd)/Data/$stage $remote_host:
EOM
done | more
