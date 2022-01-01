#!/usr/bin/env Rscript

############################################################################
#   CNC-EMDiff differential peak analysis
#   Jason W. Bacon
#   Spring, 2020
#
#   Based on OpticRegen by Andrea Rau
#   https://github.com/andreamrau/OpticRegen_2019
#   Thanks to Paul Auer for extensive guidance
############################################################################

# Choose CCA for chondro or NCA for neuro
args <- commandArgs(trailingOnly=TRUE)
if ( (length(args) != 1) || ((args[1] != "CCA") && (args[1] != "NCA")) )
{
    print("Usage: RScript diff-anal.R CCA|NCA")
    stop()
}
cell_type <- args[1]

# Compensate for misnamed sequence files.
swap_replicate_and_condition <- TRUE

# Use DiffBind summits option to recenter and trim peaks?
diffbind_summits<-FALSE

# Loading libraries can take 20 seconds or more
library(DiffBind)
library(dplyr)
library(DESeq2)

# Pause until user presses enter
pause <- function()
{
    # Change to TRUE for interactive run
    if ( FALSE )
    {
	cat("Press Enter to continue...")
	invisible(b <- scan("stdin", character(), nlines=1, quiet=TRUE))
    }
}

# Prevent print() from spewing huge tables to stdout
options(max.print=60)

# Make all output group-writable
Sys.umask(mode = 007)

setwd("10-diff-anal")

Peaks <- paste0("../9-process-peaks/p10-",
		cell_type, "-501-merged.bed")
if ( ! file.exists(Peaks) ) {
    print(paste0("Error: File ", Peaks))
    print("does not exist.  Run 9-process-peaks.sbatch first.")
    stop()
}

print("Peaks")
print(Peaks)

## Read in peaksets, identify consensus
# %>% = pipe
# Example CNC-EMDiff filename: CCA1A_S1_L001-nodup-uniq.bam
# "1A": 1 = replicate, A = time point
# Due to data mislabling, 1 = time point (condition) and A is replicate
# SampleID <- strsplit(dir("ALIGNED_TRANS/"), split=".", fixed=TRUE) %>%
# Could this just as well use the merged BAMs in 8-?
# Was using 4-bwa-mem.  Because deduped BAMs were previously stored there?
SampleID <- strsplit(dir("../6-remove-duplicates/",
    pattern=paste0(cell_type, ".*-nodup-uniq.bam")),
    split="-", fixed=TRUE) %>% lapply(., function(x) x[1]) %>%
    unlist() %>% unique()

print("SampleID")
print(SampleID)

if ( ! exists("SampleID") ) {
    print("Error: Could not generate SampleID from ../4-bwa-mem/*.bam")
    stop()
}

if ( swap_replicate_and_condition )
{
    Condition <- substr(SampleID, 4, 4) %>% unlist()    # Flip condition/replicate
} else {
    Condition <- substr(SampleID, 5, 5) %>% unlist()
}

print("Condition")
print(Condition)

if ( swap_replicate_and_condition )
{
    Replicate <- substr(SampleID, 5, 5)%>% unlist()     # Flip condition/replicate
} else {
    Replicate <- substr(SampleID, 4, 4)%>% unlist()
}

print("Replicate")
print(Replicate)
pause()

# Select only .bam files
bamReads <- paste0("../4-bwa-mem/",
		   dir("../4-bwa-mem/", pattern=paste0(cell_type, ".*.bam$")))

print("bamReads")
print(bamReads)

# rep = repeat
PeakCaller <- rep("narrow", length(SampleID))
print("PeakCaller")
print(PeakCaller)
pause()

# Create new table from existing vectors
# Specific header names required by dba.count()
samples <- data.frame(SampleID, Condition, Replicate, bamReads,
		      Peaks, PeakCaller)
print("samples")
print(samples)
pause()

readcounts_filename <- paste0(cell_type, "-readcounts_pvalsort.RData")
if ( file.exists(readcounts_filename) ) {
    print("Using saved dba.count() results...")
    print(paste("Remove ", readcounts_filename,
	" before running this script if anything has changed"))
    load(readcounts_filename)
} else {
    print("Running dba(sampleSheet=samples)")
    start_time <- Sys.time()
    peaksets_pvalsort <- dba(sampleSheet=samples)
    end_time <- Sys.time()
    print("dba() time:")
    print(end_time - start_time)
    
    # Debug
    pdf(paste0("peaksets-pvalsort-overall-clustering-", cell_type, ".pdf"))
    plot(peaksets_pvalsort)
    dev.off()   # close() pdf
    
    ## Count overlapping reads (did not recenter around peaks)
    # Count bam pileups over peaks
    # dba.count() takes a long time so use saved results if available
    # Remove the file before running this script if anything has changed
    print("Running dba.count(peaksets_pvalsort)")
    start_time <- Sys.time()
    # FIXME: Should we be using the "summits" option here?
    if ( diffbind_summits ) {
	readcounts_pvalsort <- dba.count(peaksets_pvalsort, summits=250,
				     bParallel=FALSE)
    } else {
	readcounts_pvalsort <- dba.count(peaksets_pvalsort, bParallel=FALSE)
    }
    end_time <- Sys.time()
    print("dba.count() time:")
    print(end_time - start_time)
    
    # debug
    pdf(paste0("readcounts-pvalsort-overall-clustering-", cell_type, ".pdf"))
    plot(readcounts_pvalsort)
    dev.off()
    
    # Save this one object from the session
    save(list="readcounts_pvalsort", file=readcounts_filename)
    #save(list=c("samples", "readcounts_pvalsort"), file="sampledata.RData")
}

print("readcounts_pvalsort:")
typeof(readcounts_pvalsort)
length(readcounts_pvalsort)
print(readcounts_pvalsort)
pause()

# Does not work: "Cannot coerc" error
# write.table(readcounts_pvalsort, file="readcounts_pvalsort.table",
#            quote=FALSE, sep="\t")

## Differential analysis (DESeq2) comparing all time points to 0
# colData <- samples is sufficient here
colData <- data.frame(samples)
print("colData:")
print(colData)
pause()

# Extract data we really need
# [[]] means list
rowData_pvalsort <- readcounts_pvalsort$peaks[[1]][,1:3]
rownames(rowData_pvalsort) <- paste0(# seq.int(nrow(rowData_pvalsort)), "-",
				     "chr", rowData_pvalsort$Chr, "-",
				     rowData_pvalsort$Start, "-",
				     rowData_pvalsort$End )
print("rowData_pvalsort:")
print(rowData_pvalsort)
pause()

# Get read counts for all peaks
# Organize as a table with each sample as a column
counts_pvalsort <- lapply(readcounts_pvalsort$peaks, function(x) x$Reads) %>%
  do.call("cbind", .)

# Recode X variable as a numeric indicator variable
# factor level variable instead of character for regression
colData$time_factor <- factor(colData$Condition)
print("colData$time_factor:")
print(colData$time_factor)
pause()

print("GRanges(rowData_pvalsort):")
print(GRanges(rowData_pvalsort))
pause()

# Organize input
dds_pvalsort <- DESeqDataSetFromMatrix(countData = counts_pvalsort,
				       colData = colData,
				       rowRanges = GRanges(rowData_pvalsort),
				       design = ~ time_factor)

# pauer 2020-03-16
if ( TRUE )
{
    library(RColorBrewer)
    library(pheatmap)

    vsd <- vst(dds_pvalsort, blind=FALSE)
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$SampleID)
    colnames(sampleDistMatrix) <- paste(vsd$SampleID)
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pdf(paste0("heatmap-Sample-Similarity-", cell_type, ".pdf"),
	height=8, width=11)
    pheatmap(sampleDistMatrix,
	     clustering_distance_rows=sampleDists,
	     clustering_distance_cols=sampleDists,
	     col=colors)
    dev.off()
}

# Run differential analysis
print("Running DESeq...")
dds_pvalsort <- DESeq(dds_pvalsort)
print("dds_pvalsort:")
print(dds_pvalsort)
pause()

# Generate data for Excel file

# According to file naming conventions.  Does not work: Almost no differential
# peaks.  Are we comparing replicates rather than time points due to a sample
# mixup?
if ( swap_replicate_and_condition )
{
    # Flip condition/replicate
    res_pvalsort_T1vsT0 <- results(dds_pvalsort, contrast = c("time_factor", "2", "1"))
    res_pvalsort_T2vsT0 <- results(dds_pvalsort, contrast = c("time_factor", "3", "1"))
    res_pvalsort_T2vsT1 <- results(dds_pvalsort, contrast = c("time_factor", "3", "2"))
} else {
    res_pvalsort_T1vsT0 <- results(dds_pvalsort, contrast = c("time_factor", "B", "A"))
    res_pvalsort_T2vsT0 <- results(dds_pvalsort, contrast = c("time_factor", "C", "A"))
    res_pvalsort_T2vsT1 <- results(dds_pvalsort, contrast = c("time_factor", "C", "B"))
}

print(paste0(cell_type, " T1 vs T0"))
summary(res_pvalsort_T1vsT0, alpha=0.05) ## Previously 104
write.table(res_pvalsort_T1vsT0, paste0(cell_type, "-T1-vs-T0.tsv"),
	    col.names = TRUE)

print(paste0(cell_type, " T2 vs T0"))
summary(res_pvalsort_T2vsT0, alpha=0.05) ## Previously 0
write.table(res_pvalsort_T2vsT0, paste0(cell_type, "-T2-vs-T0.tsv"),
	    col.names = TRUE)

print(paste0(cell_type, " T2 vs T1"))
summary(res_pvalsort_T2vsT1, alpha=0.05)
write.table(res_pvalsort_T2vsT1, paste0(cell_type, "-T2-vs-T1.tsv"),
	    col.names = TRUE)
