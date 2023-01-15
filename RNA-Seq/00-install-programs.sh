#!/bin/sh -e

##########################################################################
#   Run before all other scripts on supported platforms.
#   Must be run by a systems manager.
##########################################################################

if which cluster-run; then
    cluster_run=cluster-run
    srun="srun --ntasks=1 --mem=1g"
    node_spec=compute
else
    cluster_run='sh -c'
    srun=''
fi

case $(uname) in
FreeBSD)
    # Install ports on all compute nodes
    printf "Root "
    su -m root -c "$cluster_run 'pkg install -y rna-seq' $node_spec"
    ;;

*)
    cat << EOM

$0: $(uname) is not yet supported.

Please consider updating $0 to support $(uname).

EOM
    exit 1
    ;;

esac

