#############################################################################
#   Description:
#  
#   Arguments:
#
#   Returns:
#
#   History: 
#   Date        Name        Modification
#   2020-09-15  Jason Bacon Begin
#############################################################################

BEGIN {
    OFS="\t";
}
($1 !~ "baseMean") {
    
    # Read same line from other two time intervals
    do
    {
	if ( getline line20 < file20 != 1 )
	{
	    printf("Error reading %s.\n", file20);
	    exit 1;
	}
	split(line20, a20);
    }   while ( a20[1] != $1 );
    
    do
    {
	if ( getline line21 < file21 != 1 )
	{
	    printf("Error reading %s.\n", file21);
	    exit 1;
	}
	split(line21, a21);
    }   while ( a21[1] != $1 );
    
    # Throw out peaks with APV = NA in any time interval
    if ( ($7 != "NA") && (a20[7] != "NA") && (a21[7] != "NA") )
    {
	name=$1;
	lfc10=$3;
	apv10=$7;
	gsub("\"", "", name);
	split(name, a, "-");
	chr=a[1];
	start=a[2];
	end=a[3];
	gsub("chr", "", chr);

	# DA is not currently output, but save code for reference
	if ( apv10 < 0.05 )
	{
	    if ( lfc10 >= 0 )
	    {
		da="up";
	    }
	    else
	    {
		da="down";
	    }
	}
	else
	{
	    da="stable";
	}
    
	lfc20=a20[3];
	apv20=a20[7];
	lfc21=a21[3];
	apv21=a21[7];
    
	# Print floats with %s to avoid formatting in intermediate files
	printf("%s\t%u\t%u\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	       chr, start, end, name,
	       lfc10, apv10, lfc20, apv20, lfc21, apv21);
    }
}
