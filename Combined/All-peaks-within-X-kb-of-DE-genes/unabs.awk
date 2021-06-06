BEGIN {
    OFS=FS;
}
$12 ~ "-$" {
    # Move the - sign back to the beginning of the number
    gsub("-", "", $12);
    $12 = "-" $12;
}
{
    print $0;
}

