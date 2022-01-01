#!/bin/sh -e

if which cluster-run; then
    cluster_run=cluster-run
    node_spec=compute
else
    cluster_run='sh -c'
fi

case $(uname) in
FreeBSD)
    py_prefix=py38
    # Install ports on all compute nodes
    printf "Root "
    su -m root -c "$cluster_run 'pkg install -y $py_prefix-cutadapt fastqc \
	$py_prefix-multiqc samtools gffread bwa biolibc-tools \
	kallisto R hisat2' $node_spec"
	;;

*)
    printf "$0: $(uname) is not yet supported.\n"
    exit 1
    ;;

esac

# R-cran-sleuth R-cran-deseq2
printf "\nInstalling R packages...\n"
srun --ntasks=1 --mem=1g Utils/install-R-packages.R
