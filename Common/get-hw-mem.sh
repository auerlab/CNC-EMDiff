#!/bin/sh -e

# FIXME: parallel --number-of-cores reports 2 on unixdev1, should be 16
# getconf NPROCESSORS_ONLN does not work on Alma8, _NPROCESSORS_ONLN does
# _NPROCESSORS_ONLN does not work on NetBSD9
# Both forms work on FreeBSD and macOS
case $(uname) in
FreeBSD)
    sysctl -n hw.realmem
    ;;

Darwin)
    sysctl -n hw.memsize
    ;;
    
*)
    printf "Unsupported OS: $(uname).\n"
    printf "Please consider adding a case for $(uname).\n"
    exit 1
    ;;

esac
