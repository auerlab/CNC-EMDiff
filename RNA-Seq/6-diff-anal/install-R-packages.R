#!/usr/bin/env Rscript

## Default repo
local({r <- getOption("repos")
       r["CRAN"] <- "https://cran.mtu.edu/"
       options(repos=r)
})

#.libPaths(c("/usr/home/bacon/R/amd64-portbld-freebsd12.0-library/3.6",.libPaths()))
#.libPaths()

install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
