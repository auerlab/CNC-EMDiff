$12 ~ "-$" {
    # Move the - sign back to the beginning of the number
    gsub("-", "", $12);
    $12 = "-" $12;
}
{
    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	$1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12);
}

