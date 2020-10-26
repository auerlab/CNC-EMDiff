#############################################################################
#   Description:
#       Filter out non-DA peaks
#
#   History: 
#   Date        Name        Modification
#   2020-09-15  Jason Bacon Begin
#############################################################################

BEGIN {
    OFS=FS;
}
($1 != "Chr") {
    apv10=$6;
    apv20=$8;
    apv21=$10;

    # DA is not currently output, but save code for reference
    if ( (apv10 < 0.05) || (apv20 < 0.05) || (apv21 < 0.05) )
    {
	print $0;
    }
}
