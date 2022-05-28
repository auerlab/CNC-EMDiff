#!/usr/bin/env Rscript

## Default repo
local({r <- getOption("repos")
       r["CRAN"] <- "https://mirror.las.iastate.edu/CRAN/"
       options(repos=r)
})

# R_LIBS_USER should be set to ~/R/library in ~/.Renviron.

# This works, but dir is not accessible by this process
# Run the script again without removing the dir and it can install packages
# dir.create("~/R/library", recursive=TRUE)

install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Pause this thought...
# Some bugs fixed and a promise of more maintenance to come, Mar 2022.
# Sleuth has not been maintained since 2019.  Use DESeq2 instead.
# Instructions at https://pachterlab.github.io/sleuth/download do not work
# install.packages("devtools")
# devtools::install_github("pachterlab/sleuth")
# https://github.com/pachterlab/sleuth/issues/259
# remotes::install_github("pachterlab/sleuth#260")
# This seems to work
# source("http://bioconductor.org/biocLite.R")
BiocManager::install("rhdf5")
BiocManager::install("pachterlab/sleuth")
BiocManager::install("biomaRt")

# Alternative for comparison
BiocManager::install("DESeq2")
