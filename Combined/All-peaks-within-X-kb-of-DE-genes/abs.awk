BEGIN {
    OFS="\t";
}
$12 < 0 {
    # Move the - sign to the end of the number
    # sort -n tolerates this and will sort by absolute value as a result
    $12 = -$12 "-";
}
{
    print $0;
}

