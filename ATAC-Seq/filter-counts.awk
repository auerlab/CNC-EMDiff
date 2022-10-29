#############################################################################
#   Description:
#  
#   Arguments:
#
#   Returns:
#
#   History: 
#   Date        Name        Modification
#   2022-06-01  Jason Bacon Begin
#############################################################################

BEGIN {
    OFS="\t"
}
{
    t10_file="Data/15-diff-anal/" cell_type "-T1-vs-T0.tsv"
    getline t10 < t10_file
    
    t20_file="Data/15-diff-anal/" cell_type "-T2-vs-T0.tsv"
    getline t20 < t20_file
    
    t21_file="Data/15-diff-anal/" cell_type "-T2-vs-T1.tsv"
    getline t21 < t21_file
    
    if ( $1 !~ "CHR" )  # Discard column names
    {
	split(t10, a10, "[ \t]")
	split(t20, a20, "[ \t]")
	split(t21, a21, "[ \t]")
	
	#########################################################
	# Sanity check: Make sure peak names match line-for-line
	
	peak_name="chr" $1 "-" $2 "-" $3  # chr-start-end
	
	# Quotes are part of the input string in the DA files
	gsub("\"", "", a10[1])
	gsub("\"", "", a20[1])
	gsub("\"", "", a21[1])
	if ( (peak_name != a10[1]) || (peak_name != a20[1]) ||
	     (peak_name != a21[1]) ) {
	    print("Error: Files out of sync.")
	    exit(1)
	}
	
	######################################################################
	# Output counts only if significant change over at least one interval
	if ( (a10[6] <= 0.05) || (a20[6] <= 0.05) || (a21[6] <= 0.05) ) {
	    print($0, a10[6], a20[6], a21[6])
	}
    }
}
