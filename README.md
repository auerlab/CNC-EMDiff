# CNC-EMDiff
Cranial Neural Crest Ectodermal Mesenchymal Differentiation

This analysis pipeline was created for the PhD thesis of Dr. Maria Replogl
in the lab of Dr. Ava Udvadia at the University of Wisconsin -- Milwaukee.

Authors: Jason W. Bacon, Dr. Paul W. L. Auer (Sleuth analysis script)

The analysis compares gene expression data from RNA-Seq and chromatin
accessibility data from ATAC-Seq at three time points in differentiating
Mus Musculus neural and chondrocyte cells in culture.  It includes standard
differential expression and differential accessibility analyses as well as
custom scripts and newly developed C code (peak-classifier) to aid in
correlating gene expression and chromatin accessibility.

Much of this analysis is modeled on the work of Dr. Andrea Rau:

https://github.com/andreamrau/OpticRegen_2019

It was originally developed in fall 2019 and spring 2020, but has since been
overhauled to simplify the scripts, make them more readable, and make them
easily adaptable to other data sets.  Adapting the scripts for a reanalysis
of the Optic Regen data with the latest tools has proven to be quite simple.

This pipeline was developed and run on a FreeBSD cluster, with some testing
also performed on a CentOS 7 cluster, both running the SLURM batch system.
It is designed to be highly portable and easily adaptable to new data sets.

The entire pipeline can be completed in a day or two using ~64 cores.

All sbatch scripts artificially set default SLURM environment variables so
that they can be tested on a single sample outside the SLURM environment.

Programs used in the scripts are assumed to be in the user's PATH.  If not,
add appropriate steps to the sbatch scripts.  ( Altering PATH, activating
environments or containers, etc. in .bashrc or other startup scripts is not
recommended, as attempting to support various programs within a single
environment is complicated and prone to cause problems. )

The install-programs scripts will install all necessary software on a FreeBSD
cluster or workstation.

See the following links for FreeBSD setup:

http://acadix.biz/desktop-installer.php
http://acadix.biz/spcm.php
