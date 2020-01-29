#!/bin/sh -e

##########################################################################
#   Script description:
#       Trim adapters and low quality ends from Lumina reads
#       Based on work of Dr. Andrea Rau:
#       https://github.com/andreamrau/OpticRegen_2019
#
#   History:
#   Date        Name        Modification
#   2019-09-11  Jason Bacon Begin
##########################################################################

##########################################################################
#   Main
##########################################################################

# Trimgalore uses 3 processes:
#   1. trim_galore (perl)
#   2. cutadapt (python)
#   3. gzip
# All must be on the same node and should be given their own core.
# ntasks < 3 will make them share a core where CPU affinity is enabled and
# will cause you to overstep your allocation otherwise.

#SBATCH --ntasks=3 --nodes=1 --mem=256

# CentOS pkgsrc
# source /etc/bashrc
# module load /sharedapps/pkg-2017Q3/etc/modulefiles/pkgsrc/2017Q3

# In TRIMMED, we have the trimmed reads (paired end mode)
# In TRIMMED-trim1, we trim an additional 1bp from the 3' end to avoid
# alignment problems with Bowtie1
mkdir -p TRIMMED TRIMMED-trim1
for r1 in ../Raw/190822_AHFN3KDRXX-Lane2/*R1*.fastq.gz; do
    r2=$(echo ${r1} | sed -e 's|R1|R2|g')
    base_name=$(basename $r1)
    stem=${base_name%%_R*.fastq.gz}
    printf "$r1 $r2 $stem\n"
    trim_galore --stringency 3 -q 20 --paired -o TRIMMED/ \
	$r1 $r2 > TRIMMED/${stem}-trash.out

    trim_galore --stringency 3 -q 20 --paired --trim1 -o TRIMMED-trim1/ \
	$r1 $r2 > TRIMMED-trim1/${s}-trash.out
done
touch trimgalore-done
