#!/bin/sh -e

if which cluster-run; then
    cluster_run=cluster-run
    srun=srun
    node_spec=compute
else
    cluster_run='sh -c'
    srun=''
fi

case $(uname) in
FreeBSD)
    py_prefix=py38
    # Install ports on all compute nodes
    # DiffBind deps: librsvg2 pkgconf
    printf "Root "
    su -m root -c "$cluster_run 'pkg install -y $py_prefix-cutadapt fastqc \
	$py_prefix-multiqc bwa samtools $py_prefix-macs2 bedtools R \
	librsvg2 pkgconf' $node_spec"
    ;;

*)
    printf "$0: $(uname) is not yet supported.\n"
    exit 1
    ;;

esac

printf "\nInstalling R packages...\n"
$srun Utils/install-R-packages.R
