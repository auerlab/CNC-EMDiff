#!/bin/sh -e

##########################################################################
#   Script description:
#       Use BWA to align reads to a genome reference.
#
#   Dependencies:
#       Requires trimmed reads and a reference genome.  Run after
#       *-trim.sh and *-reference.sh.
##########################################################################

set -x
jobs=2

# Tried GNU parallel and ran into bugs.  Xargs just works.
echo chondro neuro \
    | xargs -n 1 -P $jobs ../../Common/redirect.sh Xargs/14-macs-peaks.sh
