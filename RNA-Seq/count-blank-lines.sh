#!/bin/sh

for file in TRIMMED/*.gz; do
    echo $file
    zcat $file | grep '^ *$' # | wc -l
done
