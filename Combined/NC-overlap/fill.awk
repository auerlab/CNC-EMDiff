#############################################################################
#   Description:
#  
#   Arguments:
#
#   Returns:
#
#   History: 
#   Date        Name        Modification
#   2021-06-02  Jason Bacon Begin
#############################################################################

BEGIN {
    OFS=",";
}
{
    for (c = 1; c <= NF; ++c) {
	if ( $c == "" ) {
	    $c = ".";
	}
    }
    print $0;
}

