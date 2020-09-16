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

$1 !~ "baseMean" && $7 != "NA" {
    name=$1;
    lfc=$3;
    apv=$7;
    gsub("\"", "", name);
    split(name, a, "-");
    chr=a[1];
    start=a[2];
    end=a[3];
    gsub("chr", "", chr);
    # DA is not currently output, but leave code just in case
    if ( apv < 0.05 )
    {
	if ( lfc >= 0 )
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
    
    # Read same line from other two time intervals
    file="../../ATAC-Seq/10-diff-anal/" cell_type "-T2-vs-T0.tsv"
    do
    {
	if ( getline t2_t0_line < file != 1 )
	{
	    printf("Error reading %s.\n", file);
	    exit 1;
	}
	split(t2_t0_line, b);
    }   while ( (b[1] ~ "baseMean") || (b[7] == "NA") );
    gsub("\"", "", b[1]);
    if ( b[1] != name )
    {
	printf("Read from %s out of sync: %s %s\n", file, name, b[1]);
	exit 1
    }
    lfc_t2_t0=b[3];
    apv_t2_t0=b[7];

    # Print floats with %s to avoid formatting in intermediate files
    printf("%s\t%u\t%u\t%s\t%s\t%s\t%s\t%s\n",
	   chr, start, end, name, lfc, apv, lfc_t2_t0, apv_t2_t0);
}
