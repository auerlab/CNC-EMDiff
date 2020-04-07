#!/bin/sh -e

##########################################################################
#   Script description:
#       Install R packages for ATAC-Seq differential analysis
#       
#   History:
#   Date        Name        Modification
#   2020-04-07  Jason Bacon Begin
##########################################################################

usage()
{
    cat << EOM

Usage: $0 prefix cran-URL

Example: $0 /raid-01/bacon/R https://cran.mtu.edu/

EOM
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 2 ]; then
    usage
fi

prefix="$1"
url="$2"

env_file=$HOME/.Renviron
profile=$HOME/.Rprofile

if ! fgrep -q R_LIBS_USER $env_file; then
    cat << EOM >> $env_file
R_LIBS_USER=
EOM
else
    printf "R_LIBS_USER already set.\n"
fi

if ! fgrep -q 'r["CRAN"]' $profile; then
    cat << EOM >> $profile
local({
  r <- getOption("repos")
  r[\"CRAN\"] <- \"$url\"
  options(repos = r)
})
EOM
else
    printf 'r["CRAN"] already set.\n'
fi

Rscript install-R-packages.R
