# CNC-EMDiff
Mus musculus neuron/chondrocyte differentiation analysis

Much of this analysis is based on the work of Andrea Rau:

https://github.com/andreamrau/OpticRegen_2019

This pipeline was developed and run on a FreeBSD cluster, with some testing
also performed on a CentOS 7 cluster.  Both clusters were running the
SLURM batch system.

The entire pipeline can be completed in a day or two using ~64 cores.

All sbatch scripts set default SLURM env variables so that they can be
tested outside SLURM environment.

Programs used in the scripts are assumed to be in the user's PATH.  If not,
add appropriate steps to the sbatch scripts.  ( Altering PATH, activating
environments or containers, etc. in .bashrc or other startup scripts is not
recommended, as attemping to support various programs within a single
environment is complicated and prone to cause problems. )

The install-programs scripts will install all necessary software on a FreeBSD
cluster or workstation.

See the following links for FreeBSD setup:

http://acadix.biz/desktop-installer.php
http://acadix.biz/spcm.php
