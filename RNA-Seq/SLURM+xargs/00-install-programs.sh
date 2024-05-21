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
    # Check for pkgsrc installed via auto-pkgsrc-setup
    if which sbatch; then
	cat << EOM

$0: You appear to be using a non-FreeBSD cluster.

You can use pkgsrc and install the biology/rna-seq package on all
compute nodes or in a shared location that compute nodes can access.

EOM
	exit 1
    else
	if which auto-pkgsrc-prefix; then
	    cd $(auto-pkgsrc-dir)/biology/rna-seq
	    bmake deinstall clean clean-depends install
	else
	    cat << EOM

$0: No pkgsrc installation found.

If you have a pkgsrc tree installed, install sysutils/auto-admin so
that $0 can use auto-pkgsrc-prefix to find it.

Otherwise, please consider updating $0 to support $(uname).

EOM
	    exit 1
	fi
    fi
    ;;

esac

