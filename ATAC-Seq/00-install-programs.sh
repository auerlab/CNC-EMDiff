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
    printf "Root "
    su -m root -c "$cluster_run 'pkg install -y atac-seq' $node_spec"
    ;;

*)
    cat << EOM

$0: $(uname) is not yet supported.

Please consider updating $0 to support $(uname).

EOM
    exit 1
    ;;

esac

