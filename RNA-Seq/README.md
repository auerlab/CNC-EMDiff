# Using the scripts

These scripts are designed to run on a SLURM HPC cluster, where all
required software is in the user's PATH.

Ensuring the availability of software may differ on each individual
HPC cluster, so it is mostly left to the user.  On FreeBSD clusters,
all software can be installed on compute nodes using
./00-install-programs.sh.  Command will be in the default PATH
on all compute nodes, so there is nothing more to do.

On many clusters, users may need to load environment modules,
conda environments, containers, etc. in order to access specific
software.  Talk to your HPC support staff for details.

## .sh vs .sbatch

Scripts with filenames ending in .sh can be run directly on the
head node, e.g. ./01-organize.sh, ./03-multiqc-raw.sh.  In some
cases, they will use SLURM's srun command to run things on a
compute node in the foreground.

Scripts with names ending in .sbatch should be submitted using the
`sbatch` command, e.g.

```
sbatch 02-qc-raw.sbatch
```

## Dependencies

Scripts can be run one at a time in numeric order for simplicity, i.e.
./00-install-programs.sh, ./01-organize.sh, etc.

Not all scripts depend on lower-numbered scripts, however, so you can
speed things up by submitting multiple scripts at the same time.
E.g. you can run 02-qc-raw.sbatch and 04-trim.sbatch at the same time,
but only after running 01-organize.sh.

Dependencies are documented in the scripts.
