#!/bin/sh -e

# Note centromere and telomere regions

printf "CDNA...\n"
cdna_ids='cdna-gene-ids.txt'
zcat ../../RNA-Seq/Reference/Mus_musculus.GRCm38.cdna.all.fa.gz \
    | awk '$4 ~ "gene:" {
	tokens = split($4, gene, ":");
	for (c = 1; c <= tokens; ++c)
	{
	    #printf("gene[%d] = %s %s\n", c, gene[c], substr(gene[c],1,4));
	    if ( substr(gene[c],1,4) == "gene" )
	    {
		split(gene[c+1], gene_id, ".");
		printf("%s\n", gene_id[1]);
		break;
	    }
	}
    }' | sort | uniq > $cdna_ids
head $cdna_ids

printf "GTF...\n"
gtf_ids='gtf-gene-ids.txt'
awk '$3 == "gene" {
    gsub("\"", "", $10);
    gsub(";", "", $10);
    printf("%s\n", $10)
    }' \
    ../../RNA-Seq/Reference/Mus_musculus.GRCm38.98.gtf \
    | sort | uniq > $gtf_ids
head $gtf_ids

printf "Diff...\n"
wc -l $gtf_ids
wc -l $cdna_ids
diff $gtf_ids $cdna_ids > ids.diff || true
egrep '^>' ids.diff > cdna-only.txt
egrep '^<' ids.diff > gtf-only.txt
wc -l cdna-only.txt
wc -l gtf-only.txt
more cdna-only.txt
more gtf-only.txt

