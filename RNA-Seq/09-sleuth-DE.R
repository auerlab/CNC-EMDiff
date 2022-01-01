#!/usr/bin/env Rscript

##########################################################################
#   Script description:
#       Sleuth analysis for CNC-EMDiff
#       
#   History:
#   Date        Name        Modification
#   Fall 2019   Paul L Auer Begin
#   2020-06-15  Jason Bacon Integrate into latest RNA-Seq pipeline
##########################################################################

############################
#############################

#######################################################
####### Our own filtering function could be ###########
####### as follows. Not using this stuff  #############
#auer_filter <- function(row, min_reads = 10, min_prop = 1)
#{
#    mean(row >= min_reads) >= min_prop
#}
#sleuth_object <- sleuth_prep(sample_to_covariates = expdesign, full_model=full_design_factor, target_mapping=t2g, filter_fun = auer_filter) 
#######################################################
######################################################

#### Organization of data, samples, replicates, etc. 
##### A = day0
####  B = day4 for Chondrocytes, day2 for Nueral Crest
####  C = day14 for Chond, day 6 for Nueral Crest

#1 - CE1A - Cd0r1
#2 - CE1B - Cd4r1
#3 - CE1C - Cd14r1
#4 - CE2A - Cd0r2
#5 - CE2B - Cd4r2
#6 - CE2C - Cd14r2
#7 - CE3A - Cd0r3
#8 - CE3B - Cd4r3
#9 - CE3C - Cd14r3
#10 - NE1A - Nd0r1
#11 - NE1B - Nd2r1
#12 - NE1C - Nd6r1
#13 - NE2A - Nd0r2
#14 - NE2B - Nd2r2
#15 - NE2C - Nd6r2
#16 - NE3A - Nd0r3
#17 - NE3B - Nd2r3
#18 - NE3C - Nd6r3


#### Within Chondro treat as a single experiment
# Read in data
# From several parameter variations early on:
# base_dir <- "4-kallisto-quant-m30-u15/"

print("Sleuth is dead.  Use DESeq2.");
quit();

# Extract numeric sample IDs from kallisto output directory names
library(stringr)
base_dir <- "Data/07-kallisto-quant";
sample_id <- str_extract(dir(base_dir),"[0-9]+");
sample_id

kal_dirs <- file.path(base_dir,dir(base_dir))
kal_dirs

expdesign <- data.frame(sample=sample_id,
			time=factor(c("T1", "T2", "T3", "T1", "T2", "T3", "T1", "T2", "T3")),
			replicate=c(1,1,1,2,2,2,3,3,3),
			path=kal_dirs,
			stringsAsFactors=FALSE)

full_design_factor <- model.matrix(formula(~ expdesign$time))
colnames(full_design_factor) <- c("CEday0", "CEday4", "CEday14")

##########################
###biomaRt stuff #########
### get gene names etc ###

# FIXME: build and release??
# Get this from GFF?
library(biomaRt)
library(dplyr)
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
	dataset = "mmusculus_gene_ensembl",
	host = 'ensembl.org');
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

library(sleuth);
so <- sleuth_prep(sample_to_covariates = expdesign, full_model=full_design_factor, target_mapping=t2g)
so <- sleuth_fit(so)
quit()


### CE Day 4 vs CE Day 0
so2 <- sleuth_wt(so, "CEday4")
res <- sleuth_results(so2, "CEday4", test_type="wt")
keep <- res[which(res$qval < 0.05),]
ce.4vs0 <- keep

### CE Day 14 vs CE Day 0
so2 <- sleuth_wt(so, "CEday14")
res <- sleuth_results(so2, "CEday14", test_type="wt")
keep <- res[which(res$qval < 0.05),]
ce.14vs0 <- keep

#### Get TPMs for each gene-transcript. Merge into results
sleuth_matrix <- data.frame(sleuth_to_matrix(so, 'obs_norm', 'tpm'))
colnames(sleuth_matrix) <- c("CE1A", "CE1B", "CE1C", "CE2A", "CE2B", "CE2C", "CE3A", "CE3B", "CE3C")
sleuth_matrix$target_id <- rownames(sleuth_matrix)

ce.4vs0 <- merge(ce.4vs0, sleuth_matrix)
ce.4vs0 <- ce.4vs0[, c(1,2,3,5,6, 14:22)]

ce.14vs0 <- merge(ce.14vs0, sleuth_matrix)
ce.14vs0 <- ce.14vs0[, c(1,2,3,5,6, 14:22)]

# mkdir() is in the R docs, but does not exist
dir.create("Sleuth-Prelim")
write.table(ce.14vs0, "Sleuth-Prelim/ce14vsce0.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(ce.4vs0, "Sleuth-Prelim/ce4vsce0.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

# Repeats for other sample groups below

## CE Day 14 vs CE Day 4
library(sleuth)
library(biomaRt)
library(dplyr)

##########################
###biomaRt stuff #########
### get gene names etc ###
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "mmusculus_gene_ensembl",
  host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
############################

# base_dir <- "~/Data/CNC-EMDiff/RNASeq/4-kallisto-quant-m30-u15/"
sample_id <- as.character(sort(as.numeric(dir(file.path(base_dir))[1:18])))[1:9]

kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
kal_dirs

expdesign <- data.frame(sample=sample_id,
			time=factor(c("A", "B", "C", "A", "B", "C", "A", "B", "C")),
			replicate=c(1,1,1,2,2,2,3,3,3),
			path=kal_dirs,
			stringsAsFactors=FALSE)

expdesign$time <- relevel(expdesign$time, ref="B")
full_design_factor <- model.matrix(formula(~ expdesign$time))
colnames(full_design_factor) <- c("CEday4", "CEday0", "CEday14")

so <- sleuth_prep(sample_to_covariates = expdesign, full_model=full_design_factor, target_mapping=t2g)
so <- sleuth_fit(so)
so2 <- sleuth_wt(so, "CEday14")
res <- sleuth_results(so2, "CEday14", test_type="wt")
keep <- res[which(res$qval < 0.05),]
ce.14vs4 <- keep

#### Get TPMs for each gene-transcript. Merge into results
sleuth_matrix <- data.frame(sleuth_to_matrix(so, 'obs_norm', 'tpm'))
colnames(sleuth_matrix) <- c("CE1A", "CE1B", "CE1C", "CE2A", "CE2B", "CE2C", "CE3A", "CE3B", "CE3C")
sleuth_matrix$target_id <- rownames(sleuth_matrix)

ce.14vs4 <- merge(ce.14vs4, sleuth_matrix)
ce.14vs4 <- ce.14vs4[, c(1,2,3,5,6, 14:22)]

write.table(ce.14vs4, "Sleuth-Prelim/ce14vsce4.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)




### Within Nuero treat as a single experiment
library(sleuth)
library(biomaRt)
library(dplyr)

##########################
###biomaRt stuff #########
### get gene names etc ###
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "mmusculus_gene_ensembl",
  host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
############################
#############################


## NE Day 2 vs NE Day 0
# base_dir <- "~/Data/CNC-EMDiff/RNASeq/4-kallisto-quant-m30-u15/"
sample_id <- as.character(sort(as.numeric(dir(file.path(base_dir))[1:18])))[10:18]

kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
kal_dirs

expdesign <- data.frame(sample=sample_id,
			time=factor(c("A", "B", "C", "A", "B", "C", "A", "B", "C")),
			replicate=c(1,1,1,2,2,2,3,3,3),
			path=kal_dirs,
			stringsAsFactors=FALSE)

expdesign$time <- relevel(expdesign$time, ref="A")
full_design_factor <- model.matrix(formula(~ expdesign$time))
colnames(full_design_factor) <- c("NEday0", "NEday2", "NEday6")


so <- sleuth_prep(sample_to_covariates = expdesign, full_model=full_design_factor, target_mapping=t2g)
so <- sleuth_fit(so)
so2 <- sleuth_wt(so, "NEday2")
res <- sleuth_results(so2, "NEday2", test_type="wt")
keep <- res[which(res$qval < 0.05),]
ne.2vs0 <- keep

## NE Day 6 vs NE Day 0
so2 <- sleuth_wt(so, "NEday6")
res <- sleuth_results(so2, "NEday6", test_type="wt")
keep <- res[which(res$qval < 0.05),]
ne.6vs0 <- keep

#### Get TPMs for each gene-transcript. Merge into results
sleuth_matrix <- data.frame(sleuth_to_matrix(so, 'obs_norm', 'tpm'))
colnames(sleuth_matrix) <- c("NE1A", "NE1B", "NE1C", "NE2A", "NE2B", "NE2C", "NE3A", "NE3B", "NE3C")
sleuth_matrix$target_id <- rownames(sleuth_matrix)

ne.2vs0 <- merge(ne.2vs0, sleuth_matrix)
ne.2vs0 <- ne.2vs0[, c(1,2,3,5,6, 14:22)]

ne.6vs0 <- merge(ne.6vs0, sleuth_matrix)
ne.6vs0 <- ne.6vs0[, c(1,2,3,5,6, 14:22)]

write.table(ne.6vs0, "Sleuth-Prelim/ne6vsne0.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
write.table(ne.2vs0, "Sleuth-Prelim/ne2vsne0.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)



### NE Day 6 vs NE Day 2
library(sleuth)
library(biomaRt)
library(dplyr)

##########################
###biomaRt stuff #########
### get gene names etc ###
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "mmusculus_gene_ensembl",
  host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
############################
#############################

# base_dir <- "~/Data/CNC-EMDiff/RNASeq/4-kallisto-quant-m30-u15/"
sample_id <- as.character(sort(as.numeric(dir(file.path(base_dir))[1:18])))[10:18]

kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
kal_dirs
quit()

expdesign <- data.frame(sample=sample_id,
			time=factor(c("A", "B", "C", "A", "B", "C", "A", "B", "C")),
			replicate=c(1,1,1,2,2,2,3,3,3),
			path=kal_dirs,
			stringsAsFactors=FALSE)

expdesign$time <- relevel(expdesign$time, ref="B")
full_design_factor <- model.matrix(formula(~ expdesign$time))
colnames(full_design_factor) <- c("NEday2", "NEday0", "NEday6")

so <- sleuth_prep(sample_to_covariates = expdesign, full_model=full_design_factor, target_mapping=t2g)
so <- sleuth_fit(so)
so2 <- sleuth_wt(so, "NEday6")
res <- sleuth_results(so2, "NEday6", test_type="wt")
keep <- res[which(res$qval < 0.05),]
ne.6vs2 <- keep

#### Get TPMs for each gene-transcript. Merge into results
sleuth_matrix <- data.frame(sleuth_to_matrix(so, 'obs_norm', 'tpm'))
colnames(sleuth_matrix) <- c("NE1A", "NE1B", "NE1C", "NE2A", "NE2B", "NE2C", "NE3A", "NE3B", "NE3C")
sleuth_matrix$target_id <- rownames(sleuth_matrix)

ne.6vs2 <- merge(ne.6vs2, sleuth_matrix)
ne.6vs2 <- ne.6vs2[, c(1,2,3,5,6, 14:22)]

write.table(ne.6vs2, "Sleuth-Prelim/ne6vsne2.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)



### For Chondro - Neuro comparisons read in the whole datasets
library(sleuth)
library(biomaRt)
library(dplyr)

##########################
###biomaRt stuff #########
### get gene names etc ###
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "mmusculus_gene_ensembl",
  host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
############################
#############################

# base_dir <- "~/Data/CNC-EMDiff/RNASeq/4-kallisto-quant-m30-u15/"
sample_id <- as.character(sort(as.numeric(dir(file.path(base_dir))[1:18])))

kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
kal_dirs

expdesign <- data.frame(sample=sample_id,
			time=factor(c("A", "B", "C", "A", "B", "C", "A", "B", "C", "D", "E", "F", "D", "E", "F", "D", "E", "F")),
			replicate=c(1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3),
			path=kal_dirs,
			stringsAsFactors=FALSE)

expdesign$time <- relevel(expdesign$time, ref="A")
full_design_factor <- model.matrix(formula(~ expdesign$time))
colnames(full_design_factor) <- c("CEday0", "CEday4", "CEday14", "NEday0", "NEday2", "NEday6")

### NE Day 0 vs CE Day 0
so <- sleuth_prep(sample_to_covariates = expdesign, full_model=full_design_factor, target_mapping=t2g)
so <- sleuth_fit(so)
so2 <- sleuth_wt(so, "NEday0")
res <- sleuth_results(so2, "NEday0", test_type="wt")
keep <- res[which(res$qval < 0.05),]
ne0.vs.ce0 <- keep

#### Get TPMs for each gene-transcript. Merge into results
sleuth_matrix <- data.frame(sleuth_to_matrix(so, 'obs_norm', 'tpm'))

colnames(sleuth_matrix) <- c("CE1A", "CE1B", "CE1C", "CE2A", "CE2B", "CE2C", "CE3A", "CE3B", "CE3C",
"NE1A", "NE1B", "NE1C", "NE2A", "NE2B", "NE2C", "NE3A", "NE3B", "NE3C")
sleuth_matrix$target_id <- rownames(sleuth_matrix)

ne0.vs.ce0 <- merge(ne0.vs.ce0, sleuth_matrix)
ne0.vs.ce0 <- ne0.vs.ce0[, c(1,2,3,5,6, 14:31)]

#### Write out results
write.table(ne0.vs.ce0, "Sleuth-Prelim/ne0vsce0.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)



##################################################################3
##### CE Day 4 vs NE Day 2
library(sleuth)
library(biomaRt)
library(dplyr)

##########################
###biomaRt stuff #########
### get gene names etc ###
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "mmusculus_gene_ensembl",
  host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
############################
#############################

# base_dir <- "~/Data/CNC-EMDiff/RNASeq/4-kallisto-quant-m30-u15/"
sample_id <- as.character(sort(as.numeric(dir(file.path(base_dir))[1:18])))

kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
kal_dirs

expdesign <- data.frame(sample=sample_id,
			time=factor(c("A", "B", "C", "A", "B", "C", "A", "B", "C", "D", "E", "F", "D", "E", "F", "D", "E", "F")),
			replicate=c(1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3),
			path=kal_dirs,
			stringsAsFactors=FALSE)

expdesign$time <- relevel(expdesign$time, ref="B")
full_design_factor <- model.matrix(formula(~ expdesign$time))
colnames(full_design_factor) <- c("CEday4", "CEday0", "CEday14", "NEday0", "NEday2", "NEday6")

so <- sleuth_prep(sample_to_covariates = expdesign, full_model=full_design_factor, target_mapping=t2g)
so <- sleuth_fit(so)
so2 <- sleuth_wt(so, "NEday2")
res <- sleuth_results(so2, "NEday2", test_type="wt")
keep <- res[which(res$qval < 0.05),]
ne2.vs.ce4 <- keep

#### Get TPMs for each gene-transcript. Merge into results
sleuth_matrix <- data.frame(sleuth_to_matrix(so, 'obs_norm', 'tpm'))

colnames(sleuth_matrix) <- c("CE1A", "CE1B", "CE1C", "CE2A", "CE2B", "CE2C", "CE3A", "CE3B", "CE3C",
"NE1A", "NE1B", "NE1C", "NE2A", "NE2B", "NE2C", "NE3A", "NE3B", "NE3C")
sleuth_matrix$target_id <- rownames(sleuth_matrix)

ne2.vs.ce4 <- merge(ne2.vs.ce4, sleuth_matrix)
ne2.vs.ce4 <- ne2.vs.ce4[, c(1,2,3,5,6, 14:31)]

#### Write out resultsv
write.table(ne2.vs.ce4, "Sleuth-Prelim/ne2vsce4.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)




######################################################
##### Finish this last part on Monday ################
######################################################


## NE Day 6 vs CE Day 14
##### Start new R-session for every new reference group below
library(sleuth)
library(biomaRt)
library(dplyr)

##########################
###biomaRt stuff #########
### get gene names etc ###
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "mmusculus_gene_ensembl",
  host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
############################
#############################

# base_dir <- "~/Data/CNC-EMDiff/RNASeq/4-kallisto-quant-m30-u15/"
sample_id <- as.character(sort(as.numeric(dir(file.path(base_dir))[1:18])))

kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
kal_dirs

expdesign <- data.frame(sample=sample_id,
			time=factor(c("A", "B", "C", "A", "B", "C", "A", "B", "C", "D", "E", "F", "D", "E", "F", "D", "E", "F")),
			replicate=c(1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3),
			path=kal_dirs,
			stringsAsFactors=FALSE)

expdesign$time <- relevel(expdesign$time, ref="C")
full_design_factor <- model.matrix(formula(~ expdesign$time))
colnames(full_design_factor) <- c("CEday14", "CEday0", "CEday4", "NEday0", "NEday2", "NEday6")

so <- sleuth_prep(sample_to_covariates = expdesign, full_model=full_design_factor, target_mapping=t2g)
so <- sleuth_fit(so)
so2 <- sleuth_wt(so, "NEday6")
res <- sleuth_results(so2, "NEday6", test_type="wt")
keep <- res[which(res$qval < 0.05),]
ne6.vs.ce14 <- keep

#### Get TPMs for each gene-transcript. Merge into results
sleuth_matrix <- data.frame(sleuth_to_matrix(so, 'obs_norm', 'tpm'))

colnames(sleuth_matrix) <- c("CE1A", "CE1B", "CE1C", "CE2A", "CE2B", "CE2C", "CE3A", "CE3B", "CE3C",
"NE1A", "NE1B", "NE1C", "NE2A", "NE2B", "NE2C", "NE3A", "NE3B", "NE3C")
sleuth_matrix$target_id <- rownames(sleuth_matrix)

ne6.vs.ce14 <- merge(ne6.vs.ce14, sleuth_matrix)
ne6.vs.ce14 <- ne6.vs.ce14[, c(1,2,3,5,6, 14:31)]

write.table(ne6.vs.ce14, "Sleuth-Prelim/ne6vsce14.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)


####### Read in the ce0vsne0, ce4vsne2, and ce14vsne6. Remove genes in ce0vsne0, and report the non-overlap for ce4vsne2 and ce14vsne6.
ce0vsne0 <- read.table("Sleuth-Prelim/ne0vsce0.txt", header=TRUE, sep="\t")
ne2vsce4 <- read.table("Sleuth-Prelim/ne2vsce4.txt", header=TRUE, sep="\t")
ne6vsce14 <- read.table("Sleuth-Prelim/ne6vsce14.txt", header=TRUE, sep="\t")

ne6vsce14.v2 <- subset(ne6vsce14, !ext_gene %in% ce0vsne0$ext_gene)
ne2vsce4.v2 <- subset(ne2vsce4, !ext_gene %in% ce0vsne0$ext_gene)

write.table(ne6vsce14.v2, "Sleuth-Prelim/ne6vsce14_v2.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(ne2vsce4.v2, "Sleuth-Prelim/ne2vsce4_v2.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

