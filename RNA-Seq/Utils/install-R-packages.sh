#!/bin/sh -e

# Not surprisingly, using dir.create() in the R script doesn't work.
# The directory is not accessible by the same Rscript process, though
# is is by subsequent runs.  :-/
mkdir -p ~/R/library
Utils/install-R-packages.R
