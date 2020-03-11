library(DiffBind)
library(dplyr)
library(DESeq2)
library(GenomicRanges)

pause <- function()
{
    cat("Press Enter to continue...")
    invisible(b <- scan("stdin", character(), nlines=1, quiet=TRUE))
}

# Update Jan 26, 2018: p-value sort on transgene alignment --------------------
# Peaks re-aligned to include transgene as well
# Use only high-quality macs2 peaks
# for(time in c(0,2,4,7,12)) {
#   macsBed <- read.delim(paste0("PEAKS_TRANS/", time, 
#                                "ATAC.nodup.unique.macs_peaks.narrowPeak"), 
#                         header=FALSE)
#   ## Sort by p-value ($V8), only keep those with p-value < 10^-10
#   keep <- which(macsBed$V8 > -log(10^-10))
#   macsBed_keep <- macsBed[keep,]
#   cat(nrow(macsBed), "...", nrow(macsBed_keep), "\n")
#   write.table(macsBed_keep, paste0("PEAKS_TRANS/", time, 
#                                    "ATAC.nodup.unique.macs_peaks.pvalsort.narrowPeak.bed"), 
#               col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
#   if(time == 0) macs_all <- macsBed_keep
#   if(time != 0) macs_all <- rbind(macs_all, macsBed_keep)
# }
# # 62703 ... 16447 
# # 46864 ... 14589 
# # 102711 ... 22400 
# # 58181 ... 17087 
# # 70846 ... 16684 
#
# ## And write out merged 
# macs_all_GR <- GRanges(seqnames=as.vector(macs_all[,1]),
#                        IRanges(start=as.numeric(as.vector(macs_all[,2])),
#                                end=as.numeric(as.vector(macs_all[,3]))),
#                        strand=rep("*", nrow(macs_all)))
# elementMetadata(macs_all_GR) <- macs_all[,-c(1:3)]
# colnames(elementMetadata(macs_all_GR)) <- c("Peak_ID","Score","Strand", 
#                                             "Fold-change", "log10pvalue", "log10qvalue",
#                                             "Relative_summit")
## Consensus peaks: 24203 peaks
# macs_merge <- reduce(macs_all_GR)
# macsGR_merge_bed <- data.frame(chr=seqnames(macs_merge),
#                                starts=start(macs_merge),
#                                end=end(macs_merge),
#                                names=c(rep(".", length(macs_merge))),
#                                scores=c(rep(".", length(macs_merge))),
#                                strands=strand(macs_merge))
# write.table(macsGR_merge_bed, paste0("PEAKS_TRANS/", 
#                                      "ALLMERGED_ATAC.nodup.unique.macs_peaks.pvalsort.narrowPeak.bed"), 
#             col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
#
# ## There is a single peak on the transgene after merging: 163:5092
# #1 99999 179 5042 0ATAC.nodup.unique.macs_peak_1 15225  . 3.18297 1527.576 1522.527 2621
# #1 99999 286 5015 2ATAC.nodup.unique.macs_peak_1 11420  . 3.70956 1147.315 1142.016 2514
# #1 99999 163 5092 4ATAC.nodup.unique.macs_peak_1 22038  . 3.30861 2208.807 2203.81 2637
# #1 99999 263 5015 7ATAC.nodup.unique.macs_peak_1 14011  . 3.04765 1406.3 1401.195 2537
# #1 99999 473 4918 12ATAC.nodup.unique.macs_peak_1 15566  . 4.09869 1562.086 1556.653 2327
# for(time in c(0,2,4,7,12)) {
#   tmp <- read.table(paste0("PEAKS_TRANS/", time, 
#                            "ATAC.nodup.unique.macs_peaks.pvalsort.narrowPeak.bed"))
#   tmp <- tmp[which(tmp$V1 == 99999),]
#   print(tmp)
# }

# Update Sep 18, 2018: peaks called on merged time pts, include pealets -------
# Col 10 is distance from peak start to summit
#macsBed <- read.delim(paste0("PEAKS_TRANS_PEAKLETS_ALLTIMES/", 
#                             "ATAC.nodup.unique.macs.peaklets_peaks.narrowPeak"), 
#                      header=FALSE)

macsBed <- read.delim("7-macs-peaklets/ATAC-CCA-macs.peaklets_peaks.narrowPeak", header=FALSE)
## Sort by p-value ($V8), only keep those with p-value < 10^-10 [44429]
# keep - index vector
# Done using awk

# Why p value and not q value?
keep <- which(macsBed$V8 > -log(10^-10))
#print("keep (rows where p < 10^-10)")
#print(keep)

# new table, keep rows in "keep"
macs_all <- macsBed[keep,]
# print("macs_all (peaks with p < 10^-10)")
# print(macs_all)

# debug: How many rows were filtered out?
cat(nrow(macsBed), "...", nrow(macs_all), "\n")

## And write out merged : AR modified to use 501 bp regions
# like bedtools - make peaklet around summit
# col 2 = start, col 10 = distance from start to peak
macs_all_GR <- GRanges(seqnames=as.vector(macs_all[,1]),
		       IRanges(start=as.numeric(as.vector(macs_all[,2])) + 
				 as.numeric(as.vector(macs_all[,10])) - 250,
			       end=as.numeric(as.vector(macs_all[,2])) + 
				 as.numeric(as.vector(macs_all[,10])) + 250),
		       strand=rep("*", nrow(macs_all)))
print("max_all_GR: 501 nt generated peaks")
print(macs_all_GR)
# line <- readline(prompt="Press [enter] to continue")
pause()

# add metadata to table
# Remove columns 1-3
elementMetadata(macs_all_GR) <- macs_all[,-c(1:3)]

# Add names for each col
colnames(elementMetadata(macs_all_GR)) <- c("Peak_ID","Score","Strand", 
					    "Fold-change", "log10pvalue",
					    "log10qvalue",
					    "Relative_summit")

## Note: we have a very small number (~50) of overlapping peaklets, I have just left
## them as is without merging (so we may be double counting a small number of reads)

# merge overlapping peaks
# Done with bedtools merge
macsGR_bed <- data.frame(chr=seqnames(macs_all_GR),
			       starts=start(macs_all_GR),
			       end=end(macs_all_GR),
			       names=c(rep(".", length(macs_all_GR))),
			       scores=c(rep(".", length(macs_all_GR))),
			       strands=strand(macs_all_GR))

## Note: some of the peaks need to be recentered (negative start)
index <- which(macsGR_bed$starts < 1)
macsGR_bed[index, "end"] <- macsGR_bed[index, "end"] - macsGR_bed[index, "starts"] + 1
macsGR_bed[index, "starts"] <- 1
## Note: some of the peaks need to be recentered (too long)
index <- which(macsGR_bed$chr == "KN149934.1" & macsGR_bed$end == 12555)
macsGR_bed[index, "end"] <- 12451; macsGR_bed[index, "starts"] <- macsGR_bed[index, "end"] - 499
index <- which(macsGR_bed$chr == "KN149962.1" & macsGR_bed$end == 136196)
macsGR_bed[index, "end"] <- 136181; macsGR_bed[index, "starts"] <- macsGR_bed[index, "end"] - 499
index <- which(macsGR_bed$chr == "KN150281.1" & macsGR_bed$end == 3781)
macsGR_bed[index, "end"] <- 3763; macsGR_bed[index, "starts"] <- macsGR_bed[index, "end"] - 499
index <- which(macsGR_bed$chr == "KN150391.1" & macsGR_bed$end == 26616)
macsGR_bed[index, "end"] <- 26543; macsGR_bed[index, "starts"] <- macsGR_bed[index, "end"] - 499
#write.table(macsGR_bed, paste0("PEAKS_TRANS_PEAKLETS_ALLTIMES/",
# "ALLMERGED_ATAC.nodup.unique.macs.peaklets_peaks.pvalsort.narrowPeak_501bp.bed")
write.table(macsGR_bed, "CCA-merged-peaks.bed",
	    col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

#print("macsGR_bed: Merged peaks")
#print(macsGR_bed)
#pause()

#------------------------------------------------------------------------------------

## Read in peaksets, identify consensus
# %>% = pipe
# https://github.com/andreamrau/OpticRegen_2019/blob/master/1_preprocessing/ATACseq/7_remove_MTduplicates.sh
# Example OpticRegen filename: ALIGNED_TRANS/12ATAC-3.bam
# Condition = 12, replicate = 3
# Example CNC-EMDiff filename: CCA1A_S1_L001-nodup-uniq.bam
# "1A": 1 = condition, A = time point
# SampleID <- strsplit(dir("ALIGNED_TRANS/"), split=".", fixed=TRUE) %>%
SampleID <- strsplit(dir("4-bwa-mem/", pattern="CCA.*.bam"), split=".", fixed=TRUE) %>%
  lapply(., function(x) x[1]) %>%
  unlist() %>%
  unique()
# 12ATAC-3

print("SampleID")
print(SampleID)

# Condition <- strsplit(SampleID, split="ATAC", fixed=TRUE) %>%
Condition <- substr(SampleID, 5, 5) %>% unlist()
# 12

print("Condition")
print(Condition)
pause()

# Replicate <- strsplit(SampleID, split="-", fixed=TRUE) %>%
Replicate <- substr(SampleID, 4, 4)%>% unlist()
# %>% as.numeric()
# 3

print("Replicate")
print(Replicate)
pause()

# Filter out bam.*, presumably leaving only .bam files
# bamReads <- paste0("ALIGNED_TRANS/",
#                dir("ALIGNED_TRANS/")[-grep("bam.", dir("ALIGNED_TRANS/"))])
bamReads <- paste0("4-bwa-mem/",
		 dir("4-bwa-mem/", pattern="CCA.*.bam$"))

print("bamReads")
print(bamReads)
pause()

# Peaks <- paste0("PEAKS_TRANS_PEAKLETS_ALLTIMES/", Condition, 
#                 "ATAC.nodup.unique.macs_peaks.pvalsort.narrowPeak.bed")
# Peaks <- paste0("PEAKS_TRANS_PEAKLETS_ALLTIMES/",
#                 "ALLMERGED_ATAC.nodup.unique.macs.peaklets_peaks.pvalsort.narrowPeak_501bp.bed")
Peaks <- "CCA-merged-peaks.bed"
print("Peaks")
print(Peaks)
pause()

# rep = repeat
PeakCaller <- rep("narrow", length(SampleID))
print("PeakCaller")
print(PeakCaller)
pause()

# Create new table from existing vectors
samples <- data.frame(SampleID, Condition, Replicate, bamReads, Peaks, PeakCaller)
print("samples")
print(samples)
pause()

readcounts_file <- "readcounts_pvalsort.RData"
if ( file.exists(readcounts_file) ) {
    print("Using saved dba.count() results...")
    print(paste("Remove ", readcounts_file,
	" before running this script if anything has changed"))
    load(readcounts_file)
} else {
    print("Running dba(sampleSheet=samples)")
    start_time <- Sys.time()
    peaksets_pvalsort <- dba(sampleSheet=samples)
    end_time <- Sys.time()
    print("dba() time:")
    print(end_time - start_time)
    pdf("peaksets_pvalsort_overall-clustering.pdf")
    # Debug?
    plot(peaksets_pvalsort)
    dev.off()   # close() pdf
    
    ## Count overlapping reads (did not recenter around peaks)
    # Count bam pileups over peaks
    # dba.count() takes a long time so use saved results if available
    # Remove the file before running this script if anything has changed

    print("Running dba.count(peaksets_pvalsort)")
    start_time <- Sys.time()
    readcounts_pvalsort <- dba.count(peaksets_pvalsort)
    end_time <- Sys.time()
    print("dba.count() time:")
    print(end_time - start_time)
    
    # debug
    pdf("readcounts_pvalsort_overall-clustering.pdf")
    plot(readcounts_pvalsort)
    dev.off()
    
    # Save this one object from the session
    # load("readcounts_pvalsort.RData") to view
    save(list="readcounts_pvalsort", file=readcounts_file)
    #save(list=c("samples", "readcounts_pvalsort"), file="sampledata.RData")
}
print("readcounts_pvalsort:")
print(readcounts_pvalsort)
pause()

# Does not work: "Cannot coerc" error
#write.table(readcounts_pvalsort, file="readcounts_pvalsort.table",
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
print("rowData_pvalsort:")
print(rowData_pvalsort)
pause()

# Get read counts for all peaks
# Organize as a table with each sample as a column
counts_pvalsort <- lapply(readcounts_pvalsort$peaks, function(x) x$Reads) %>%
  do.call("cbind", .)
#print("counts_pvalsort:")
#print(counts_pvalsort)

# Recode X variable as a numeric indicator variable
# factor level variable instead of character for regression
colData$time_factor <- factor(colData$Condition)
print("colData$time_factor:")
print(colData$time_factor)
pause()

# Organize input
dds_pvalsort <- DESeqDataSetFromMatrix(countData = counts_pvalsort,
				       colData = colData,
				       rowRanges = GRanges(rowData_pvalsort),
				       design = ~ time_factor)

# Run differential analysis
print("Running DESeq...")
dds_pvalsort <- DESeq(dds_pvalsort)
print("dds_pvalsort:")
print(dds_pvalsort)
pause()

# Adjust this for Maria's time points:
# 
# Generate data for Excel file
## 184 total are differential (down from 208 in previous analysis)
res_pvalsort_BvsA <- results(dds_pvalsort, contrast = c("time_factor", B, A))
res_pvalsort_CvsA <- results(dds_pvalsort, contrast = c("time_factor", C, A))
res_pvalsort_CvsB <- results(dds_pvalsort, contrast = c("time_factor", C, B))

summary(res_pvalsort_BvsA, alpha=0.05) ## Previously 104
# out of 42198 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 141, 0.33% 
# LFC < 0 (down)   : 9, 0.021% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 14726, 35% 
# (mean count < 31)
summary(res_pvalsort_CvsA, alpha=0.05) ## Previously 0
# out of 42198 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 1, 0.0024% 
# LFC < 0 (down)   : 0, 0% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 0, 0% 
# (mean count < 1)
summary(res_pvalsort_CvsB, alpha=0.05)
# out of 42198 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 0, 0% 
# LFC < 0 (down)   : 0, 0% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 0, 0% 
# (mean count < 1)
stop()

# Create CSVs with differntially accessible peaks
## No differential peaks at 7 to 0
for(i in c(2,4,12)) {
  res_pvalsort <- get(paste0("res_pvalsort_", i, "vs0"))
  sigpeaks_pvalsort <- data.frame(rowRanges(dds_pvalsort[which(get(paste0("res_pvalsort_", i, "vs0"))$padj < 0.05),])[,1:2],
				    get(paste0("res_pvalsort_", i, "vs0"))[which(get(paste0("res_pvalsort_", i, "vs0"))$padj < 0.05),])
  o <- order(abs(sigpeaks_pvalsort$log2FoldChange), decreasing=TRUE)
  sigpeaks_pvalsort_order <- sigpeaks_pvalsort[o,]
    sigpeaks_pvalsort_order$new_seqnames <- paste0(sigpeaks_pvalsort_order$seqnames,
						 ":", sigpeaks_pvalsort_order$start,
						 "-", sigpeaks_pvalsort_order$end)
  allpeaks_pvalsort <- data.frame(rowRanges(dds_pvalsort)[,1:2],
				  get(paste0("res_pvalsort_", i, "vs0")))
  allpeaks_pvalsort$new_seqnames <- paste0(allpeaks_pvalsort$seqnames,
					   ":", allpeaks_pvalsort$start,
					   "-", allpeaks_pvalsort$end)
  saveRDS(sigpeaks_pvalsort_order, paste0("DIFFERENTIAL_vs0/sigpeaks_pvalsort_order_", i, "vs0_peaklets_501bp.rds"))
  saveRDS(allpeaks_pvalsort, paste0("DIFFERENTIAL_vs0/allpeaks_pvalsort_", i, "vs0_peaklets_501bp.rds"))
    saveRDS(res_pvalsort, paste0("DIFFERENTIAL_vs0/res_pvalsort_", i, "vs0_peaklets_501bp.rds"))
  write.csv(sigpeaks_pvalsort_order, paste0("DIFFERENTIAL_vs0/res_pvalsort_", i, "vs0_peaklets_501bp.csv"))
}


saveRDS(dds_pvalsort, "DIFFERENTIAL_vs0/dds_pvalsort_peaklets_501bp.rds")
save.image("DIFFERENTIAL_vs0/differential-analysis_diffBind_pvalsort_trans_peaklets_501bp.RData")

# peaks opening or closing from time 0
## Count number of overall differential / nondifferential
tmp <- cbind(res_pvalsort_2vs0$padj, res_pvalsort_4vs0$padj,
	     res_pvalsort_7vs0$padj, res_pvalsort_12vs0$padj)
tmp2 <- cbind(res_pvalsort_2vs0$log2FoldChange, res_pvalsort_4vs0$log2FoldChange,
	     res_pvalsort_7vs0$log2FoldChange, res_pvalsort_12vs0$log2FoldChange)
tmp <- tmp < 0.05
tmp[is.na(tmp)] <- 0
table(rowSums(tmp))

## Count number of open/closed differential peaks
colSums((tmp2 < 0) * tmp)
colSums((tmp2 > 0) * tmp)
colSums(tmp)

table(rowSums((tmp2 < 0) * tmp))
table(rowSums((tmp2 > 0) * tmp))

tmp <- cbind(res_pvalsort_2vs0$padj, res_pvalsort_4vs0$padj,
	     res_pvalsort_7vs0$padj, res_pvalsort_12vs0$padj)
tmp[is.na(tmp)] <- 1
tmp2 <- cbind(res_pvalsort_2vs0$log2FoldChange, res_pvalsort_4vs0$log2FoldChange,
	      res_pvalsort_7vs0$log2FoldChange, res_pvalsort_12vs0$log2FoldChange)
tmp[c(11499, 11719),]
tmp2[c(11499, 11719),]
rowData_pvalsort[c(11499, 11719),]


#-----------------------------------------------------
## Output bed files of differential ATAC peaks

tmp_all <- data.frame(chr=character(), start=numeric(), end=numeric(),
		      stringsAsFactors = FALSE)
for(i in c(2,4,12)) {
  tmp <- readRDS(paste0("DIFFERENTIAL_vs0/sigpeaks_pvalsort_order_", i, 
			"vs0_peaklets_501bp.rds"))
  tmp_bed <- tmp %>%
    select(chr=seqnames, start, end)
  write.table(tmp_bed, file = paste0("DIFFERENTIAL_vs0/sigpeaks_", i, 
				     "vs0_peaklets_501bp.bed"),
	      row.names=FALSE, col.names=FALSE, quote=FALSE)
  tmp_all <- bind_rows(tmp_all, tmp_bed)
}

tmp_all <- unique(tmp_all)
write.table(tmp_all, file = paste0("DIFFERENTIAL_vs0/sigpeaks_all_peaklets_501bp.bed"),
	    row.names=FALSE, col.names=FALSE, quote=FALSE)


#-----------------------------------------------------
## Output bed files of directional differential ATAC peaks


for(i in c(2,12)) {
      tmp <- readRDS(paste0("DIFFERENTIAL_vs0/sigpeaks_pvalsort_order_", i, 
			"vs0_peaklets_501bp.rds"))
  for(type in c("up", "down")) {
    tmp_bed <- tmp %>%
      mutate(sign = ifelse(log2FoldChange > 0, "up", "down")) %>%
      filter(sign == type) %>%
      select(chr=seqnames, start, end)
    print(dim(tmp_bed))
    write.table(tmp_bed, file = paste0("DIFFERENTIAL_vs0/sigpeaks_", i, 
				       "vs0_peaklets_501bp_", type, ".bed"),
		row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  }

}


tmp_2 <- readRDS(paste0("DIFFERENTIAL_vs0/sigpeaks_pvalsort_order_", 2, 
			"vs0_peaklets_501bp.rds")) %>%
  select(new_seqnames) %>% unlist()
tmp_12 <- readRDS(paste0("DIFFERENTIAL_vs0/sigpeaks_pvalsort_order_", 12, 
			"vs0_peaklets_501bp.rds")) %>%
  select(new_seqnames) %>% unlist()
tmp <- c(tmp_2, tmp_12)
tmp[duplicated(tmp)]

