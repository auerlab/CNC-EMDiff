#############################################################################
#   Description:
#  
#   Arguments:
#
#   Returns:
#
#   History: 
#   Date        Name        Modification
#   2022-11-07  Jason Bacon Begin
#############################################################################

BEGIN {
    agree=0;
    disagree=0;
}
{
    if ( $1 == "Feature" ) {
	print $0;
    }
    #printf("%s %s\n", $5, $9);
    if ( ($1 ~ "ENS") && (($5 < 0.05 && $9 >= 0.05) || ($5 >= 0.05 && $9 < 0.05)) ) {
	++disagree;
	print $0;
    }
    else {
	++agree;
    }
}
END {
    printf("Agree: %d  Disagree: %d\n", agree, disagree);
}

