#!/bin/sh -e

##########################################################################
#   Find overlap of DA peaks from neuro and chondro that are within 250
#   kb of commonly regulated DE genes
#
#     1. report % overlapping peaks
#     2. report % overlapping peaks with matching temporal accessibility
#        patterns
#     3. for each overlap report DA temporal accessibility patterns for
#        chondro and neuro  2. Find overlap of DA peaks from neuro and chondro that are 
##########################################################################

# Filter peaks for genes in genes-common*.txt and differential accessibility
# across any time interval

