#!/bin/sh -e

if which cluster-run; then
    cluster_run=cluster-run
    srun=srun
    node_spec=compute
else
    cluster_run='sh -c'
    srun=''
fi

py_prefix=py38
progs="$py_prefix-cutadapt fastqc $py_prefix-multiqc bwa samtools \
	$py_prefix-macs2 bedtools R librsvg2 pkgconf"

case $(uname) in
FreeBSD)
    # Install ports on all compute nodes
    # DiffBind deps: librsvg2 pkgconf
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

printf "\nInstalling R packages...\n"
$srun Utils/install-R-packages.R
