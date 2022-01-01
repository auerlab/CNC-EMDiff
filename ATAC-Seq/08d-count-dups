#!/bin/sh -e

##########################################################################
#   Script description:
#       Count duplicate reads removed
#       
#   History:
#   Date        Name        Modification
#   2021-02-18  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 0 ]; then
    usage
fi

cat << EOM

Below are the deduplicated percentages reported by FastQC on the raw SAM files.
These are an estimate the % of reads remaining after removing duplicates.
They are only an estimate based on sequences present in the first 100,000
reads:

https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3 Analysis Modules/8 Duplicate Sequences.html

Some additional information not clear from the FastQC docs can be found here:

https://www.biostars.org/p/107402/

EOM
cd 5-qc-sam
for file in *.zip; do
    dir=${file%.zip}
    printf "$dir\t"
    unzip -oq $file $dir/fastqc_data.txt
    fgrep '#Total Deduplicated Percentage' $dir/fastqc_data.txt | awk '{ print $4 }'
    #rm -rf $dir
done
cd ..

cat << EOM

Below is an estimate of duplicated reads based on counting BWA alignments
before and after deduplication.

EOM

printf "Raw-file\tDedup-file\tRaw-alignments\tDedup-alignments\t%%\n"
for sam in 4-bwa-mem/*.sam; do
    sample=$(basename $sam)
    sample=${sample%.sam}
    bam=$(ls 6-remove-duplicates/$sample-nodup-uniq.bam)
    raw_alignments=$(samtools view $sam | wc -l)
    dedup_alignments=$(samtools view $bam | wc -l)
    printf "%s\t%s\t%u\t%u\t%u\n" $(basename $sam) $(basename $bam) \
	$raw_alignments $dedup_alignments \
	$(($dedup_alignments * 100 / $raw_alignments))
done
