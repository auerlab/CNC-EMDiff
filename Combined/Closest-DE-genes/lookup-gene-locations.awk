#############################################################################
#   Description:
#       Find genes in a file sorted by ID and extract location info
#
#   Usage:
#       awk -f this-script -v gene_locations_file=gene-locations-file gene-list-file
#
#   Arguments:
#       Main file arg is the file containing gene names to look up
#       Variable gene_locations_file names file with gene id, chromosome,
#       and positions, generally extracted from a GTF
#
#   History: 
#   Date        Name        Modification
#   2020-06-22  Jason Bacon Begin
#############################################################################

{
    gene_id=$1;
    # Skip lines in gene-names file until we find one matching the gene_id
    while ( (getline < gene_locations_file ) && ($4 < gene_id) )
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

