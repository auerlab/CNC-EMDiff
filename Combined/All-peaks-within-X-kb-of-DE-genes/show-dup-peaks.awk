#############################################################################
#   Description:
#       
#   Arguments:
#
#   Returns:
#
#   History: 
#   Date        Name        Modification
#   2020-09-22  Jason Bacon Begin
#############################################################################

BEGIN {
    last_pos = 0;
    dup_lines = 0;
    OFS="\t";
}
{
    while ( $2 == last_pos )
    {
	if ( dup_lines++ == 0 )
	{
	    printf("\n%s %s %s\n", NR, last_pos, dup_lines);
	    print "last_line = " last_line;
	}
	print "This line = " $0;
	getline;
    }
    last_line = $0;
    last_pos = $2;
    dup_lines = 0;
}
