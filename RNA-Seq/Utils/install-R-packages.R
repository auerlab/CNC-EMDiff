#!/usr/bin/env Rscript

## Default repo
local({r <- getOption("repos")
       r["CRAN"] <- "https://mirror.las.iastate.edu/CRAN/"
       options(repos=r)
})

# R_LIBS_USER should be set in ~/.Renviron.

install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Sleuth has not been maintained since 2019.  Use DESeq2 instead.
# Instructions at https://pachterlab.github.io/sleuth/download do not work
# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")
# install.packages("devtools")
# devtools::install_github("pachterlab/sleuth")
# This seems to work
# BiocManager::install("pachterlab/sleuth")

BiocManager::install("DESeq2")
