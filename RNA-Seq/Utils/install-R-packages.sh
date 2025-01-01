#!/bin/sh -e

# R_LIBS_USER should be set to shared-path/R/library in ~/.Renviron
r_libs_user=~/R/library
printf "R_LIBS_USER=$r_libs_user\n" >> ~/.Renviron

# Not surprisingly, using dir.create() in the R script doesn't work.
# The directory is not accessible by the same Rscript process, though
# is is by subsequent runs.  :-/
mkdir -p $r_libs_user
Utils/install-R-packages.R
