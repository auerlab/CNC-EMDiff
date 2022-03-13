#!/usr/bin/env Rscript

## Default repo
local({r <- getOption("repos")
       r["CRAN"] <- "https://mirror.las.iastate.edu/CRAN/"
       options(repos=r)
})

# R_LIBS_USER should be set to ~/R/library in ~/.Renviron.

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

if (!requireNamespace("ChIPpeakAnno", quietly = TRUE))
    BiocManager::install("ChIPpeakAnno")

if (!requireNamespace("AnnotationHub", quietly = TRUE))
    BiocManager::install("AnnotationHub")
