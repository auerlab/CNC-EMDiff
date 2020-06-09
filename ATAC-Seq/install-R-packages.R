#!/usr/bin/env Rscript

if(!require(dplyr, quietly = TRUE))
    install.packages("dplyr")

if(!require(RColorBrewer, quietly = TRUE))
    install.packages("RColorBrewer")

if(!require("pheatmap", quietly = TRUE))
    install.packages("pheatmap")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")

if (!requireNamespace("GenomicRanges", quietly = TRUE))
    BiocManager::install("GenomicRanges")

if (!requireNamespace("DiffBind", quietly = TRUE))
    BiocManager::install("DiffBind")
