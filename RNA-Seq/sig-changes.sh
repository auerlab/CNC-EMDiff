#!/bin/sh -e

##########################################################################
#   Description:
#       Help the user identify biologically significant fold-changes
#       in the DE output.  This involves filtering based on
#       user-provided criteria, showing fold-change, P-value, other
#       statistics, available gene function information, etc.
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
    printf "Usage: $0 \n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 0 ]; then
    usage
fi

# Possibly use GAF files directly, or tools sets like goatools or ermineJ
# wget http://geneontology.org/gene-associations/goa_human.gaf.gz

