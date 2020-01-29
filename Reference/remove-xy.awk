#############################################################################
#   Description:
#       Remove X and Y chromosomal info from mouse cDNA fasta file
#
#   History: 
#   Date        Name        Modification
#   2019-10-05  Jason Bacon Begin
#############################################################################

{
    # As long as the description line indicates and X or Y chromosome,
    # discard it and everything to the next sequence line.
    while ( $0 ~ "GRCm38.[XYxy]" )
    {
	do
	{
	    getline
	}   while ( $0 !~ "^>" );
    }
    print $0
}

