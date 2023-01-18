#!/bin/sh -e

for cell_type in chondro neuro; do
    awk -v cell_type=$cell_type -f filter-counts.awk \
	Data/15-diff-anal/$cell_type-counts.tsv \
	> Data/15-diff-anal/$cell_type-counts-significant.tsv
done
