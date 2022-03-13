#!/bin/sh -e

# Record OS and software versions
uname -a > Logs/16-peak-classes/uname.txt

build=$(../Common/genome-build.sh)
release=$(../Common/genome-release.sh)

cat << EOM

build = $build
release = $release

EOM

# FIXME: Why was 100 used?
# gff=Mus_musculus.GRCm$build.100.gff3.gz
gff=Mus_musculus.GRCm$build.$release.gff3.gz

if [ ! -e $gff ]; then
    for fetch in fetch curl wget; do
	if which $fetch; then
	    break
	fi
    done
    if [ 0$fetch = 0 ]; then
	printf "No fetch program found (tried fetch, curl, wget).\n"
	exit 1
    fi
    if [ $fetch = curl ]; then
	fetch="$fetch -O"
    fi
    $fetch ftp://ftp.ensembl.org/pub/release-$release/gff3/mus_musculus/$gff
fi
output_dir=12-peak-classes
mkdir -p $output_dir

which peak-classifier

##########################################################################
#   All peaks
##########################################################################

for peaks_file in Data/14-process-peaks/p10-*-501-merged.bed; do
    printf "\n===\nMACS2 output:   $peaks_file\n"
    
    overlaps_file=$output_dir/$(basename ${peaks_file%.bed}-overlaps.tsv)
    printf "Overlaps file:  $overlaps_file\n"
    peak-classifier --midpoints $peaks_file $gff $overlaps_file
    
    filtered_overlaps=${overlaps_file%.tsv}-filtered.tsv
    printf "Filtered:       $filtered_overlaps\n"
    filter-overlaps $overlaps_file $filtered_overlaps five_prime_utr intron \
	exon upstream1000 upstream10000 upstream100000 upstream-beyond
done

##########################################################################
#   DA peaks
##########################################################################

pval=0.05
for file in Data/15-diff-anal/*-T*.tsv; do
    printf "\n===\nDESeq2 output:     $file\n"
    peaks_file=$output_dir/$(basename ${file%.tsv}-p$pval.bed)
    printf "Converting to bed: $peaks_file\n"
    fgrep -v baseMean $file | tr '"' ' ' | \
	awk -v pval=$pval '$7 < pval {
	    split($1, a, "-");
	    # R write.table() outputs numbers like 22000000 as 2.2e+7 even
	    # when embedded in the chromosome name.  Use %u to convert to
	    # integer format for bedtools.
	    printf("%s\t%u\t%u\t%s-%u-%u\t0\n",
		substr(a[1],4,1),a[2],a[3],a[1],a[2],a[3]);
	}' > $peaks_file

    overlaps_file=${peaks_file%.bed}-overlaps.tsv
    peak-classifier --midpoints $peaks_file $gff $overlaps_file
    
    filtered_overlaps=${overlaps_file%.tsv}-filtered.tsv
    filter-overlaps $overlaps_file $filtered_overlaps five_prime_utr intron \
	exon upstream1000 upstream10000 upstream100000 upstream-beyond
done
