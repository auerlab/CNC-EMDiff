#!/bin/sh -e

if [ $(whoami) == root ]; then
    printf "$0 should not be run as root.  It installs CRAN packages in $HOME.\n"
    exit 1
fi

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
    # Install ports on all compute nodes
    # DiffBind deps: librsvg2 pkgconf
    # coreutils: gsort is much faster than native sort and peak-classifier
    # will use it if available.
    py_prefix=py38
    progs="$py_prefix-cutadapt fastqc $py_prefix-multiqc bwa samtools \
	    $py_prefix-macs2 bedtools R librsvg2 pkgconf biolibc-tools \
	    peak-classifier coreutils"
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

install_dir=~/R/library
mkdir -p $install_dir

if ! grep '^R_LIBS_USER' ~/.Renviron; then
    cat << EOM >> ~/.Renviron
R_LIBS_USER=$install_dir
# R_LIBS=~/R/library
EOM
fi

printf "\nInstalling R packages...\n"
$srun Utils/install-R-packages.sh
