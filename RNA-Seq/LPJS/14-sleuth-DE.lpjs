#!/bin/sh -e

#
#   Dependencies:
#       Requires kallisto abundances.  Run after *-kallisto-quant.sbatch.

# FIXME: Just guessing.  Monitor and adjust
#SBATCH --output=Logs/14-sleuth-DE/slurm-%A_%a.out
#SBATCH --error=Logs/14-sleuth-DE/slurm-%A_%a.err
#SBATCH --mem=2g

./sleuth-DE.R
