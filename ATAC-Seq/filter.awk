#############################################################################
#   Description:
#       Filter out features from one bed file that are not listed in another
#
#   Usage:
#       awk -f this-script -v file_to_filter=big-file little-file
#
#   Arguments:
#       Main file arg is the file containing only the features to keep
#       Variable file_to_filter specifies the name of the file to be filtered
#
#   History: 
#   Date        Name        Modification
#   2020-02-20  Jason Bacon Begin
#############################################################################

{
    name=$4;
    # Skip features in big file until we find one matching the name in little file
    while ( (getline < file_to_filter) && ($4 != name) )
    {
    }
    # Print the matching feature
    printf("%s\t%s\t%s\t%s\t0\n", $1, $2, $3, "chr"$1"-"$2"-"$3);
}

