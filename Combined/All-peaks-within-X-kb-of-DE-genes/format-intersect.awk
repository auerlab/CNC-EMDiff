#############################################################################
#   Description:
#       Format output of bedtools intersect
#
#   History: 
#   Date        Name        Modification
#   2020-09-16  Jason Bacon Begin
#############################################################################

BEGIN {
    OFS="\t";
}
{
    # Peak info
    chr=$1;
    peak_start=$2;
    peak_end=$3;
    peak_name=$4;
    lfc10=$5;
    apv10=$6;
    lfc20=$7;
    apv20=$8;
    lfc21=$9;
    apv21=$10;
    
    # Gene info
    vicinity_start=$12;
    vicinity_end=$13;
    gene_name=$14;
    tss=$15;
    if ( tss < peak_start ) {
	gene_distance = tss - peak_start;
    }
    else if ( tss > peak_end ) {
	gene_distance = tss - peak_end;
    }
    else {
	gene_distance = 0;
    }
    # Format final floating point output
    printf("%s\t%u\t%u\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%d\n",
	   chr, peak_start, peak_end, peak_name,
	   lfc10, apv10, lfc20, apv20, lfc21, apv21,
	   gene_name, gene_distance);
}
