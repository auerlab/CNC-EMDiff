#############################################################################
#   Description:
#  
#   Arguments:
#
#   Returns:
#
#   History: 
#   Date        Name        Modification
#   2023-01-15  Jason Bacon Begin
#############################################################################

$1 == "id:" && $2 ~ "GO:" {
    id=$2
    getline
    name=$0
    sub("name: ", "", name);
    printf("%s\t%s\n", id, name);
}

