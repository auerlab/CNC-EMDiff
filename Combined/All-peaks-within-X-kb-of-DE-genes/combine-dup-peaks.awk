#############################################################################
#   Description:
#       Combine duplicate peak lines into one containing all genes and their
#       distances from the peak
#
#   History: 
#   Date        Name        Modification
#   2020-09-24  Jason Bacon Begin
#############################################################################

BEGIN {
    last_pos = 0;
    OFS=FS;
}
{
    if ( NR == 1 )
    {
	printf("%s", $0);
    }
    else if ( $2 != last_pos )
    {
	printf("\n%s", $0);
    }
    else
    {
	while ( $2 == last_pos )
	{
	    printf("\t%s\t%s", $11, $12);
	    getline;
	}
	printf("\n%s", $0);
    }
    last_pos = $2;
}
