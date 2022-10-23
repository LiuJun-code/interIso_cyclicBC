#!/bin/bash
. /home/local/CSI/jl95muce/openfoam/OpenFOAM-v2012/etc/bashrc WM_NCOMPROCS=2; export WM_COMPILE_OPTION=Debug
/usr/bin/gdb "$@"