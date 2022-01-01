#!/bin/sh -e

if which cluster-run; then
    cluster_run=cluster-run
    node_spec=compute
else
    cluster_run='sh -c'
fi

py_prefix=py38
progs="$py_prefix-cutadapt fastqc $py_prefix-multiqc samtools gffread bwa biolibc-tools kallisto R hisat2"

case $(uname) in
FreeBSD)
    # Install ports on all compute nodes
    printf "Root "
    su -m root -c "$cluster_run 'pkg install -y $progs' $node_spec"
    ;;

*)
    cat << EOM

$0: $(uname) is not yet supported.

Please find another way to install

$progs

or consider updating $0 to support $(uname).

EOM
    exit 1
    ;;

esac

# R-cran-sleuth R-cran-deseq2
printf "\nInstalling R packages...\n"
srun --ntasks=1 --mem=1g Utils/install-R-packages.R
