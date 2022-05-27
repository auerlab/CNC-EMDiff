#!/bin/sh -e

if [ $(whoami) == root ]; then
    printf "$0 should not be run as root.  It installs CRAN packages in $HOME.\n"
    exit 1
fi

if which cluster-run; then
    cluster_run=cluster-run
    srun="srun --ntasks=1 --mem=1g"
    node_spec=compute
else
    cluster_run='sh -c'
    srun=''
fi

# libgit2 is a dep for R devtools
py_prefix=py38
progs="$py_prefix-cutadapt fastq-trim fastqc $py_prefix-multiqc samtools gffread bwa biolibc-tools kallisto R hisat2 libgit2"

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

install_dir=~/R/library
mkdir -p $install_dir

if ! grep '^R_LIBS_USER' ~/.Renviron; then
    cat << EOM >> ~/.Renviron
R_LIBS_USER=$install_dir
# R_LIBS=~/R/library
EOM
fi

printf "\nInstalling R packages...\n"
$srun Utils/install-R-packages.R
