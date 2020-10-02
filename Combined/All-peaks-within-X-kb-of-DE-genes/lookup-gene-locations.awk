#############################################################################
#   Description:
#       Find genes in a file sorted by ID and extract location info
#
#   Usage:
#       awk -f this-script -v gene_locations=gene-locations-file gene-list-file
#
#   Arguments:
#       Main file arg is the file containing all gene IDs to look up
#       Variable gene_locations contains gene id, chromosome, and positions
#
#   History: 
#   Date        Name        Modification
#   2020-06-22  Jason Bacon Begin
#############################################################################

BEGIN {
    OFS=FS;
}
{
    gene_id=$1;
    # Skip lines in gene-ids file until we find one matching the gene_id
    while ( (getline < gene_locations ) && ($4 < gene_id) )
    {
    }
    if ( $4 == gene_id )
    {
	printf("%s\t%s\t%s\t%s\n", $1, $2, $3, gene_id);
    }
    else
    {
	# Indicate gene ID not found in GTF
	printf("-1\t-1\t-1\t%s\n", gene_id);
    }
}

