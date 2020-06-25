# CNC-EMDiff
Mus musculus neuron/chondrocyte differentiation analysis

Much of this analysis is based on the work of Andrea Rau:

https://github.com/andreamrau/OpticRegen_2019

This pipeline was developed and run primarily on FreeBSD 12.1 with some
testing on a CentOS 7 cluster.

Sbatch scripts are designed to run under the SLURM workload manager.  The
entire pipeline can be completed in a day or two using ~64 cores.
The wrapper script "slurm-sim" can be used to run sbatch scripts on a
workstation or laptop, which of course will take much longer.

The install-programs script will install all necessary software on a FreeBSD
cluster or workstation. See the following links for FreeBSD setup:

http://acadix.biz/desktop-installer.php
http://acadix.biz/spcm.php
