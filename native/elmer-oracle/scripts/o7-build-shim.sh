#!/usr/bin/env bash
# Build the O7 kernel shim against one pristine Elmer arm (report §5.2
# item 8). Oracle-only: the .so links that arm's libelmersolver and is
# loaded by `elmer-oracle o7-replay` and the gated `tests/o7.rs`.
#
# usage: o7-build-shim.sh <build-dir> <out.so>
set -euo pipefail
build=$(cd "$1" && pwd)
out=$2
src=$(dirname "$0")/../o7/o7_shim.f90
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
# The arm's own Fortran flags, for a uniform build. The shim only copies
# values in and out of the kernels, so its flags cannot change a result.
flags="$(sed -n 's/^CMAKE_Fortran_FLAGS:STRING=//p' "$build/CMakeCache.txt") $(sed -n 's/^CMAKE_Fortran_FLAGS_RELEASE:STRING=//p' "$build/CMakeCache.txt")"
# libelmersolver's RUNPATH names the container build tree, so its
# dependencies do not resolve from here; the Rust loader (`o7.rs`) preloads
# them from <build-dir> by absolute path before it opens the shim.
gfortran $flags -fPIC -shared -J "$work" -I "$build/fmodules" "$src" \
  -o "$out" -L "$build/fem/src" -lelmersolver
echo "built $out against $build ($flags)"
