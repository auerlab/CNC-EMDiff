#!/bin/sh -e

##########################################################################
#   Script description:
#       Locate ATAC-Seq peaks within ${kb}kb of DE gene TSS for each
#       gene cluster.
#       Proximal peaks are defined as within 1kb
#       Distal peaks are defined as within ${kb}kb
#       
#   History:
#   Date        Name        Modification
#   2020-06-12  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 KB-from-TSS\n"
    exit 1
}


##########################################################################
#   Function description:
#       Pause until user presses return
##########################################################################

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}

##########################################################################
#   Function description:
#       Print a header to an output file
#       
#   History:
#   Date        Name        Modification
#   2020-07-03  Jason Bacon Begin
##########################################################################

header()
{
    printf "Chr\tStart\tEnd\tPeak-name\tT1-T0-LFC\tT1-T0-APV\tT2-T0-LFC\tT2-T0-APV\tT2-T1-LFC\tT2-T1-APV\tGene-name\tGene-Distance\n"
    return 0
}


##########################################################################
#   Function description:
#       
#   Arguments:
#       
#   Returns:
#       
#   History:
#   Date        Name        Modification
#   2020-10-05  Jason Bacon Begin
##########################################################################

check_for_zeroes()
{
    local file=$1
    
    # Check for 0s
    $awk '$6 == 0 || $8 == 0 || $10 == 0' \
	$file > $tmpdir/zero.tsv
    if [ $(cat $tmpdir/zero.tsv | wc -l) != 0 ]; then
	printf "Error: Found zeroes in $pd_peaks_file.\n"
	printf "Check formatting statements.\n"
	exit 1
    fi
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
kb=$1

awk=$(../../Common/find-awk.sh)
results_dir=Proximal+distal-${kb}kb
mkdir -p $results_dir

bedtools --version > $results_dir/bedtools-version.txt

# Sort by chromosome and start position.  We should not have to look at end
# position since peaks should not overlap as this point.
bed_sort="sort -k1,1 -k2,2n"

# Sort by chromosome, start position, and abs(distance from TSS), with
# some help from awk filters below.
abs_sort="sort -k1,1 -k2,2n -k12,12n"

##########################################################################
#   Extract gene locations from Ensembl GTF.  They are also in the CDNA
#   reference, but CDNA does not document features like exons, UTRs, etc.
#   Ensembl GTF has many more gene IDs represented than CDNA reference
#   and CDNA contains a small number of gene IDs not in the GTF.
#   all-but-xy.transcripts.clean.fa is generated from either GTF or CDNA
#   Last run used CDNA, so there will be some gene IDs in the clusters
#   not found in the GTF.  Other parts of this analysis will need to look
#   at exons, so use GTF here for consistency.
##########################################################################

./cleanup.sh

rna_ref_dir="../../RNA-Seq/Reference"
genome_build=$(../../Common/genome-build.sh)
genome_release=$(../../Common/genome-release.sh)
gtf="$rna_ref_dir/Mus_musculus.GRCm$genome_build.$genome_release.gtf"
if [ ! -e $gtf ]; then
    (cd $(dirname $gtf) && ./fetch-gtf.sh)
fi

tmpdir=Temp
mkdir -p $tmpdir
gene_locations=$tmpdir/gene-locations-$genome_release.tsv
if [ ! -e $gene_locations ]; then
    printf "Generating $gene_locations...\n"
    # Output format: chr start stop gene-ID
    # Sort by gene ID for awk script that follows
    $awk '$1 >= 1 && $1 <= 19 && $3 == "gene" && $13 == "gene_name" {
	chr=$1;
	start=$4;
	end=$5;
	gene_name=$14;
	# Strip quotes and ; from gene name in GTF
	gsub("\"", "", gene_name);
	gsub(";", "", gene_name);
	printf("%s\t%u\t%u\t%s\n", chr, start, end, gene_name);
    }' $gtf | sort -k 4 > $gene_locations
fi
wc -l $gene_locations
head -5 $gene_locations

##########################################################################
#   Get gene IDs from RNA-Seq clustering results (provided by Maria)
#   Clustering performed by Paul Auer
##########################################################################

for cell_type in CCA NCA; do
    
    # Col 16 of augmented Sleuth output contains cluster # for each gene
    # List of unique values in this column enumerates clusters
    tsv="../../RNA-Seq/Clusters/$cell_type.tsv"
    clusters=$($awk '$1 != "target_id" { print $16 }' $tsv | sort | uniq)
    
    for cluster in $clusters; do
	printf "=== $cell_type, cluster $cluster ===\n"
	
	# Generate list of all gene names represented in this cluster
	cluster_gene_list=$tmpdir/genes-$cell_type-c$cluster.txt
	# Maria's hand-made cluster spreadsheets
	# $16 = cluster-number, $3 = gene name
	$awk "\$16 == $cluster { print \$3 }" $tsv | sort | uniq \
	    > $cluster_gene_list
	
	# Generate BED file with genes in this cluster and their locations
	cluster_gene_locations=${cluster_gene_list%.txt}.bed
	$awk -v gene_locations=$gene_locations -f lookup-gene-locations.awk \
	    $cluster_gene_list > $cluster_gene_locations

	# Generate bed file with just chr, TSS, and gene name
	# Filter out gene IDs not found in GTF (tagged with -1 by
	# lookup-gene-locations.awk) and replace end with start since we're
	# measuring distance from TSS in either direction.
	$awk '$1 != -1 {
		chr=$1;
		start=$2;
		gene_name=$4;
		printf("%s\t%s\t%s\t%s\t%s\n", chr, start, start+1, gene_name, start);
	    }' $cluster_gene_locations > $tmpdir/temp.bed

	# Generate BED file with gene TSS position +/-${kb}kb
	# Peaks that overlap this region are within ${kb}kb of the gene
	cluster_tss_vicinity=${cluster_gene_locations%.bed}-${kb}kb.bed
	bedtools slop \
	    -b ${kb}000 \
	    -g $rna_ref_dir/chromosome-sizes.tsv \
	    -i $tmpdir/temp.bed \
	    | $bed_sort > $cluster_tss_vicinity
	wc -l $cluster_tss_vicinity
	
	#####################################################################
	#   Find all peaks (both differential and stable accessibility)
	#   within ${kb}kb of DE gene TSS by cluster
	#####################################################################
	
	printf "Finding proximal and distal peaks (<= ${kb}kb from TSS)...\n"
	# Generate matching BED file from DESeq2 output
	# All intervals T1-vs-T0, T2-vs-T0, T2-vs-T1 contain the same
	# set of peaklets, so pick one from which to extract common info
	# and read LFC and APV from the other two.
	peaks_file=$tmpdir/peaks-$cell_type-c$cluster.tsv
	# DESeq2 annoyingly decides to use scientific notation for
	# boundaries sometimes, so convert boundaries to plain integers
	# while massaging to BED format.
	# Filter out header line containing "baseMean"
	da_dir=../../ATAC-Seq/10-diff-anal
	file10=$da_dir/$cell_type-T1-vs-T0.tsv
	file20=$da_dir/$cell_type-T2-vs-T0.tsv
	file21=$da_dir/$cell_type-T2-vs-T1.tsv
	$awk -f gen-peaks-file.awk -v file20=$file20 -v file21=$file21 \
	    $file10 | $bed_sort >> $peaks_file
	wc -l $peaks_file
	
	# Find peaks within intersecting TSS+/-${kb}kb
	pd_peaks_file=$results_dir/$(basename ${cluster_gene_locations%.tsv}-proximal+distal-peaks.tsv)
	header > $pd_peaks_file
	# See printfs above to verify:
	# Cols 1-10:    Peak             chr,start,end,name,LFC,APV,LFC,APV
	# Cols 11-15:   TSS-vicinity     chr,vstart,vend,gene-name,tss
	bedtools intersect -sorted -wa -wb -a $peaks_file \
	    -b $cluster_tss_vicinity | $bed_sort | uniq \
	    | $awk -f format-intersect.awk \
	    | format-cols '\t' 0.01 3 5 6 7 8 9 10 >> $pd_peaks_file
	wc -l $pd_peaks_file
	check_for_zeroes $pd_peaks_file
	
	printf "Generating variable-column file...\n"
	var_col_peaks_file=${pd_peaks_file%.tsv}-var-cols.tsv
	# Number of columns should agree with printf in format-intersect.awk
	$awk -f abs.awk $pd_peaks_file | $abs_sort | $awk -f unabs.awk \
	    | $awk -f combine-dup-peaks.awk > $var_col_peaks_file
	
	printf "Filtering for 750 nt peaks...\n"
	wide_peaks_file=${pd_peaks_file%.tsv}-over-750.tsv
	header > $wide_peaks_file
	$awk '$3 - $2 > 750 { print $0 }' $pd_peaks_file \
	    >> $wide_peaks_file
	check_for_zeroes $wide_peaks_file

	printf "Generating variable-column file for 750 nt peaks...\n"
	var_col_peaks_file=${wide_peaks_file%.tsv}-var-cols.tsv
	# Number of columns should agree with printf in format-intersect.awk
	$awk -f abs.awk $wide_peaks_file | $abs_sort | $awk -f unabs.awk \
	    | $awk -f combine-dup-peaks.awk > $var_col_peaks_file
	
	printf "Filtering for DA peaks...\n"
	da_peaks_file=${pd_peaks_file%.tsv}-da.tsv
	header > $da_peaks_file
	$awk -f da-peaks.awk $pd_peaks_file >> $da_peaks_file
	check_for_zeroes $da_peaks_file

	printf "Generating variable-column file for DA peaks...\n"
	var_col_peaks_file=${da_peaks_file%.tsv}-var-cols.tsv
	# Number of columns should agree with printf in format-intersect.awk
	$awk -f abs.awk $da_peaks_file | $abs_sort | $awk -f unabs.awk \
	    | $awk -f combine-dup-peaks.awk > $var_col_peaks_file
	
	#####################################################################
	#   Generate fasta file for each set of proximal+distal peaks
	#   seqkit rmdup should be useless here since we did sort|uniq
	#   above, but include it as a redundant check.
	#####################################################################
	
	printf "Generating fastas...\n"
	atac_ref_dir=../../ATAC-Seq/Reference
	genome_fasta=$atac_ref_dir/Mus_musculus.GRCm$genome_build.dna.autosomes.fa

	fasta_file=$results_dir/$(basename ${peaks_file%.bed}.fasta)
	$awk '$1 != "Chr"' $pd_peaks_file \
	    | bedtools getfasta -fi $genome_fasta -bed - \
	    | seqkit rmdup > $fasta_file
	wc -l $fasta_file

	fasta_file=$results_dir/$(basename ${peaks_file%.bed}-over-750.fasta)
	$awk '$1 != "Chr"' $wide_peaks_file \
	    | bedtools getfasta -fi $genome_fasta -bed - \
	    | seqkit rmdup > $fasta_file
	wc -l $fasta_file

	fasta_file=$results_dir/$(basename ${peaks_file%.bed}-da.fasta)
	$awk '$1 != "Chr"' $da_peaks_file \
	    | bedtools getfasta -fi $genome_fasta -bed - \
	    | seqkit rmdup > $fasta_file
	wc -l $fasta_file
    done
done
