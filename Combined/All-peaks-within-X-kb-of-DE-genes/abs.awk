$12 < 0 {
    # Move the - sign to the end of the number
    # sort -n tolerates this and will sort by absolute value as a result
    $12 = -$12 "-";
}
{
    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	$1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12);
}

