#!/bin/sh -e

##########################################################################
#   Script description:
#       Organize files and create directories needed for analysis
#
#       Link raw files to standardized names for that clearly describe
#       conditions and replicates, and are easily parsed by subsequent
#       scripts.
#
#       Researchers involved in sample-prep generally don't think about
#       how filename conventions impact bioinformatics analysis, so this
#       simple step can avoid confusion throughout the pipeline.  In
#       addition, linking this way can correct for sample mixups, etc.
#
#       Use links to preserve the original files and document the mapping.
#
#   Dependencies:
#       Creates required directory structure.
#       Run this before all other scripts.
#       
#   History:
#   Date        Name        Modification
#   2021-09-25  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 0 ]; then
    usage
fi

mkdir -p Results Logs
scripts=$(ls 0[2-9]-* 1[0-9]-*)
for script in $scripts; do
    stage=${script%.*}
    mkdir -p Results/$stage Logs/$stage
done

##############################################################################
# RNA-Seq:
#
# The digit in the 3rd character of the filename (1,2,3) indicates replicate
# The letter in the 4th character (A, B, C) represents time point
#   A = time 0, B = time 1, C = time 2
# The number following _S is the sample #
#
# Samples 1-9 are Chondrocyte and 10-18 are Neural
# Every third sample (those with the same letter A, B, or C in the 4th
# character of the filename) is from the same time point, e.g.
#
# CE1A_S1, CE2A_S4, and CE3A_S7 are chondrocyte time 0
# CE1B_S2, CE2B_S5, and CE3B_S8 are Chondrocyte time 1
#
# CE1A_S1_L002-R1.fastq.xz
# CE1A_S1_L002-R2.fastq.xz
# CE1B_S2_L002-R1.fastq.xz
# CE1B_S2_L002-R2.fastq.xz
# CE1C_S3_L002-R1.fastq.xz
# CE1C_S3_L002-R2.fastq.xz
# CE2A_S4_L002-R1.fastq.xz
# CE2A_S4_L002-R2.fastq.xz
# CE2B_S5_L002-R1.fastq.xz
# CE2B_S5_L002-R2.fastq.xz
# CE2C_S6_L002-R1.fastq.xz
# CE2C_S6_L002-R2.fastq.xz
# CE3A_S7_L002-R1.fastq.xz
# CE3A_S7_L002-R2.fastq.xz
# CE3B_S8_L002-R1.fastq.xz
# CE3B_S8_L002-R2.fastq.xz
# CE3C_S9_L002-R1.fastq.xz
# CE3C_S9_L002-R2.fastq.xz
# NE1A_S10_L002-R1.fastq.xz
# NE1A_S10_L002-R2.fastq.xz
# ...
# NE3C_S18_L002-R2.fastq.xz
##############################################################################

# CE1A_S1_L002-R1.fastq.xz
cd Results
rm -rf 01-organize/Raw-renamed
mkdir -p 01-organize/Raw-renamed
cd 01-organize/Raw-renamed
for rep in 1 2 3; do
    for time in A B C; do
	numeric_time=$(echo $time | tr "ABC" "123")
	for read in 1 2; do
	    orig=$(ls ../../../../../Raw/190822_AHFN3KDRXX-Lane2-RNA/CE${rep}${time}_*_L002_R${read}_001.fastq.xz)
	    ls -l $orig
	    sample=$(echo $orig | awk -F '_' '{ print $3 }' | sed -e 's|S|sample|')
	    ln -sf $orig \
		chondro-$sample-rep${rep}-time${numeric_time}-R${read}.fastq.xz
	    orig=$(ls ../../../../../Raw/190822_AHFN3KDRXX-Lane2-RNA/NE${rep}${time}_*_L002_R${read}_001.fastq.xz)
	    sample=$(echo $orig | awk -F '_' '{ print $3 }' | sed -e 's|S|sample|')
	    ln -sf $orig \
		neuro-$sample-rep${rep}-time${numeric_time}-R${read}.fastq.xz
	done
    done
done
