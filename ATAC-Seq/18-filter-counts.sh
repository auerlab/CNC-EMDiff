#!/bin/sh -e

for cell_type in chondro neuro; do
    awk -v cell_type=$cell_type -f filter-counts.awk \
	Results/15-diff-anal/$cell_type-counts.tsv \
	> Results/15-diff-anal/$cell_type-counts-significant.tsv
done
