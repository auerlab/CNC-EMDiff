#!/bin/sh -e

# FIXME: Just guessing.  Monitor and adjust
#SBATCH --output=Logs/13-sleuth-DE/slurm-%A_%a.out
#SBATCH --error=Logs/13-sleuth-DE/slurm-%A_%a.err
#SBATCH --mem=2g

./sleuth-DE.R
