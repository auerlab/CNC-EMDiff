#!/bin/sh -e

fasda=chondro-time1-time2-FC.txt
sleuth=Data/13-sleuth-DE/ch-t2-vs-t1.txt
wc $fasda $sleuth

header="$(head -1 $fasda | awk '{ printf("%20s %7s %7s %6s %5s\n", $1, $2, $3, $4, $5); }')    SC1    SC2  SFC    SPV"
printf "$header\n"
for feature in $(awk '$1 ~ "ENS" { print $1 }' $sleuth); do
    printf "$(awk -v feature=$feature '$1 == feature { printf("%20s %7.1f %7.1f %6.1f %0.3f\n", $1, $2, $3, $4, $5); }' $fasda) "
    #fasda_pval=$(awk -v feature=$feature '$1 == feature { print $5 }' $fasda)
    sleuth_pval=$(awk -v feature=$feature '$1 == feature { print $4 }' $sleuth)
    sleuth_t1=$(awk -v feature=$feature '$1 == feature { print $6 + $9 + $12 }' $sleuth)
    sleuth_t2=$(awk -v feature=$feature '$1 == feature { print $7 + $10 + $13 }' $sleuth)
    if [ $sleuth_t1 = 0 ]; then
	sleuth_fc=0.0
    else
	sleuth_fc=$(echo "$sleuth_t2 / $sleuth_t1" | bc -l)
    fi
    printf "%6.1f %6.1f %4.1f %0.4f\n" \
	$sleuth_t1 $sleuth_t2 $sleuth_fc $sleuth_pval
done
