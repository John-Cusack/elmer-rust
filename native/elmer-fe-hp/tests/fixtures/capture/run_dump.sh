#!/bin/bash
# usage: run_dump.sh <case-dir> <out.jsonl>   (serial, niced; pristine binary only)
B=/home/john/runs/elmer-fork/orig_src/build/fem
export OMP_NUM_THREADS=1 HWLOC_COMPONENTS=-gl LD_LIBRARY_PATH=$B/src ELMER_HOME=$B/src ELMER_LIB=$B/src
export DUMP_OUT=$2
[ -n "$NOEDGE" ] || export NOEDGE=
W=/home/john/firstmate/data/elmer-swarm-l15/refs
( cd "$1" && nice -n 10 gdb -batch -ex "source $W/dump.py" --args $B/src/ElmerSolver > "$2.gdb.log" 2>&1 )
