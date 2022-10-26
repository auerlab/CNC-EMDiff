#!/bin/sh -e

fasda=chondro-time1-time2-FC.txt
sleuth=Data/13-sleuth-DE/ch-t2-vs-t1.txt
wc $fasda $sleuth

for feature in $(awk '$1 ~ "ENS" { print $1 }' $sleuth); do
    printf "$(awk  -v feature=$feature '$1 == feature' $fasda) "
    #fasda_pval=$(awk -v feature=$feature '$1 == feature { print $5 }' $fasda)
    sleuth_pval=$(awk -v feature=$feature '$1 == feature { print $4 }' $sleuth)
    printf "$sleuth_pval\n"
done
