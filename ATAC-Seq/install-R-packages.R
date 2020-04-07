#!/usr/bin/env Rscript

install.packages("dplyr")
install.packages("RColorBrewer")
install.packages("pheatmap")

BiocManager::install("DESeq2")
BiocManager::install("GenomicRanges")
BiocManager::install("DiffBind")
