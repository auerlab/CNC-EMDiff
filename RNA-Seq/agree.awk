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
    min=0.03;
    max=0.07;
}
{
    if ( $1 == "Feature" ) {
	print $0;
    }
    if ( ($1 ~ "ENS") && (($5 < min && $9 >= max) || ($5 >= max && $9 < min)) ) {
	++disagree;
	print $0;
    }
    else {
	++agree;
    }
}
END {
    printf("Min = %f  Max = %f\n", min, max);
    printf("Agree = %d  Disagree = %d  Agreement = %d%%\n",
	    agree, disagree, agree * 100 / (agree + disagree));
}

