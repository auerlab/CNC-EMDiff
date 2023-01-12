#!/bin/sh -e

mkdir -p Data Logs
scripts=$(ls 0[2-9]-* [1-9][0-9]-*)
for script in $scripts; do
    stage=${script%.*}
    mkdir -p Data/$stage Logs/$stage
done

##############################################################################
# ATAC-Seq:
#
# Was meant to use the same naming convention as RNA-Seq, but replicate and
# time point were swapped in the sample naming.
#
# The digit in the 3rd character of the filename (1,2,3) indicates time point
# The letter in the 4th character (A, B, C) represents replicate
#   A = rep 0, B = rep 1, C = rep 2
# The number following _S is the sample #
#
# Samples 1-9 are Chondrocyte and 10-18 are Neural
# Every third sample (those with the same letter A, B, or C in the 4th
# character of the filename) is from the same replicate, e.g.
#
# CCA1A_S1, CCA2A_S4, and CCA3A_S7 are chondrocyte rep 0
# CCA1B_S2, CCA2B_S5, and CCA3B_S8 are chondrocyte rep 1
#
# CCA1A_S1_L001-R1.fastq.xz
# CCA1A_S1_L001-R2.fastq.xz
# CCA1B_S2_L001-R1.fastq.xz
# CCA1B_S2_L001-R2.fastq.xz
# CCA1C_S3_L001-R1.fastq.xz
# CCA1C_S3_L001-R2.fastq.xz
# CCA2A_S4_L001-R1.fastq.xz
# CCA2A_S4_L001-R2.fastq.xz
# CCA2B_S5_L001-R1.fastq.xz
# CCA2B_S5_L001-R2.fastq.xz
# CCA2C_S6_L001-R1.fastq.xz
# CCA2C_S6_L001-R2.fastq.xz
# CCA3A_S7_L001-R1.fastq.xz
# CCA3A_S7_L001-R2.fastq.xz
# CCA3B_S8_L001-R1.fastq.xz
# CCA3B_S8_L001-R2.fastq.xz
# CCA3C_S9_L001-R1.fastq.xz
# CCA3C_S9_L001-R2.fastq.xz
# NCA1A_S10_L001-R1.fastq.xz
# NCA1A_S10_L001-R2.fastq.xz
# ...
# NCA3C_S18_L001-R2.fastq.xz
##############################################################################

# CCA1A_S1_L001-R1.fastq.xz
cd Data
rm -rf 01-organize/Raw-renamed
mkdir -p 01-organize/Raw-renamed
cd 01-organize/Raw-renamed

for time in 1 2 3; do
    for rep in A B C; do
	numeric_rep=$(echo $rep | tr "ABC" "123")
	for read in 1 2; do
	    orig=$(ls ../../../../Raw/190822_AHFN3KDRXX-Lane1-ATAC/CCA${time}${rep}_*_L001_R${read}_001.fastq.xz)
	    sample=$(echo $orig | awk -F '_' '{ print $3 }' | sed -e 's|S|sample|')
	    ln -sf $orig \
		chondro-$sample-rep${numeric_rep}-time${time}-R${read}.fastq.xz
	    orig=$(ls ../../../../Raw/190822_AHFN3KDRXX-Lane1-ATAC/NCA${time}${rep}_*_L001_R${read}_001.fastq.xz)
	    sample=$(echo $orig | awk -F '_' '{ print $3 }' | sed -e 's|S|sample|')
	    ln -sf $orig \
		neuro-$sample-rep${numeric_rep}-time${time}-R${read}.fastq.xz
	done
    done
done
