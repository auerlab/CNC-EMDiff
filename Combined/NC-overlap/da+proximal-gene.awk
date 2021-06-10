#############################################################################
#   Description:
#       Filter peaks file for peaks with differential accessibility and
#       proximity to a list of genes
#
#   $6 = T1-T0 APV
#   $8 = T2-T0 APV
#   $10 = T2-T1 APV
#
#   History: 
#   Date        Name        Modification
#   2021-06-08  Jason Bacon Begin
#############################################################################

# In case we want to look at LFC
function abs(x)
{
    return x < 0 ? -x : x;
}

# If any APV is < 0.05, it's DA
($6 < 0.05) || ($8 < 0.05) || ($10 < 0.05) {
    chrom=$1;
    start=$2;
    end=$3;
    gene=$11;
    split(genestr, genes);
    # printf("gene = %s\n", gene);
    for (c in genes) {
	if ( gene == genes[c] ) {
	    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
		    chrom, start, end, gene, $6, $8, $10);
	}
    }
}
