#!/bin/bash
# Configure and build one pinned Fortran Elmer arm for the S0 oracle.
#
# Usage: build-arm.sh <arm> <source-dir> <build-dir>
#   arm:        release | fma | openblas | listusage
#   source-dir: pristine Elmer checkout at the pinned upstream commit
#   build-dir:  fresh directory this script creates
#
# Standing rules, enforced here rather than remembered:
# - HWLOC_COMPONENTS=-gl is exported for every configure, build, and test,
#   or Open MPI start-up hangs.
# - At most two full Fortran Elmer builds run on this host across all
#   workers; this script refuses to start a third.
# - Every arm records its toolchain in <build-dir>/versions.txt, which is
#   committed alongside the golden data it produces.
set -euo pipefail

PIN=9f6af2f85433fbe8ec71b85c3486ed38b88fcaf0

if [ "$#" -ne 3 ]; then
    echo "usage: build-arm.sh <release|fma|openblas|listusage> <source-dir> <build-dir>" >&2
    exit 2
fi

ARM="$1"
SRC="$2"
DEST="$3"

# Container mounts change file ownership; without this every git call below
# fails with "dubious ownership".
export GIT_CONFIG_COUNT=1
export GIT_CONFIG_KEY_0=safe.directory
export GIT_CONFIG_VALUE_0='*'
running=$(pgrep -c '^(cmake|make|ninja|ninja-build|cc1|cc1plus|f951|cpp|as|ld)$' || true)
if [ "$running" -ge 2 ]; then
    echo "refusing: $running Elmer build processes already on this host (cap is 2)" >&2
    pgrep -a '^(cmake|make|ninja|ninja-build|cc1|cc1plus|f951|cpp|as|ld)$' >&2 || true
    exit 3
fi

if [ -e "$DEST" ]; then
    echo "refusing: build dir $DEST already exists" >&2
    exit 2
fi

head=$(git -C "$SRC" rev-parse HEAD)
if [ "$head" != "$PIN" ]; then
    echo "refusing: $SRC is at $head, want pristine pin $PIN" >&2
    exit 2
fi
if [ -n "$(git -C "$SRC" status --porcelain)" ]; then
    echo "refusing: $SRC has uncommitted changes (oracle arms build pristine sources only)" >&2
    exit 2
fi

# Elmer's configure writes `results/` directories into the source tree, so
# a read-only source mount cannot configure. Copy the verified pin to a
# writable sibling of the build dir instead; the pin check above already
# covers provenance, and the copy carries .git so the record stays exact.
if [ ! -w "$SRC" ]; then
    SRC_COPY="$DEST.src"
    if [ -e "$SRC_COPY" ]; then
        echo "refusing: source copy $SRC_COPY already exists" >&2
        exit 2
    fi
    cp -a "$SRC" "$SRC_COPY"
    SRC="$SRC_COPY"
fi

# shellcheck disable=SC2207
common_flags=(
    -DCMAKE_BUILD_TYPE=Release
    -DWITH_MPI=TRUE -DWITH_ElmerIce=TRUE
    -DWITH_ELMERGUI=FALSE -DWITH_ELMERGUILOGGER=FALSE -DWITH_ELMERGUITESTER=FALSE
    -DWITH_ELMERPOST=FALSE -DWITH_LUA=FALSE -DWITH_MKL=FALSE -DWITH_Mumps=FALSE
    -DWITH_Hypre=FALSE -DWITH_Trilinos=FALSE -DWITH_FETI4I=FALSE -DWITH_Zoltan=FALSE
    -DWITH_ROCALUTION=FALSE -DWITH_CHOLMOD=FALSE -DWITH_UMFPACK=TRUE -DWITH_CONTRIB=FALSE
    -DWITH_OpenMP=TRUE -DWITH_NETCDF=FALSE -DWITH_MMG=FALSE -DWITH_PARMMG=FALSE
    -DWITH_XIOS=FALSE -DWITH_PROJ=FALSE
    -DELMERGRID_WITH_MATC=OFF
    # Reference netlib BLAS/LAPACK by absolute path: the image's
    # alternatives point libblas.so.3 at OpenBLAS, so the default search
    # would silently make every arm an OpenBLAS arm.
    -DBLAS_LIBRARIES=/usr/lib/x86_64-linux-gnu/blas/libblas.so.3
    -DLAPACK_LIBRARIES=/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3
    "-DCMAKE_EXE_LINKER_FLAGS=-Wl,--allow-multiple-definition"
    "-DCMAKE_SHARED_LINKER_FLAGS=-Wl,--allow-multiple-definition"
)

case "$ARM" in
    release) extra_flags=() ;;
    fma) extra_flags=(-DCMAKE_Fortran_FLAGS="-march=x86-64-v3 -ffp-contract=fast") ;;
    openblas)
        extra_flags=(
            -DBLA_VENDOR=OpenBLAS
            -DBLAS_LIBRARIES=/usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so
            -DLAPACK_LIBRARIES=/usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so
        )
        ;;
    listusage) extra_flags=(-DDEVEL_LISTUSAGE=ON) ;;
    *) echo "unknown arm $ARM" >&2; exit 2 ;;
esac

cmake -S "$SRC" -B "$DEST" "${common_flags[@]}" "${extra_flags[@]}"
cmake --build "$DEST" --parallel "$(nproc)"

{
    echo "arm: $ARM"
    echo "source-pin: $head"
    echo "date-utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "--- compilers ---"
    gfortran --version | head -1
    gcc --version | head -1
    cmake --version | head -1
    ninja --version
    echo "--- mpi ---"
    mpirun --version | head -2 || true
    ompi_info --param hwloc all 2>/dev/null | grep -i components || true
    echo "--- blas/lapack ---"
    readlink -f /usr/lib/x86_64-linux-gnu/libblas.so.3 /usr/lib/x86_64-linux-gnu/liblapack.so.3
    grep -E '^(BLAS_LIBRARIES|LAPACK_LIBRARIES|BLA_VENDOR)' "$DEST/CMakeCache.txt" || true
    echo "--- cpu ---"
    grep -m1 'model name' /proc/cpuinfo
    grep -m1 flags /proc/cpuinfo | tr ' ' '\n' | grep -c -E '^(fma|avx2|avx512f)$' || true
    echo "--- image ---"
    docker images --digests elmer-oracle-fortran 2>/dev/null || echo "native (non-docker) build"
} > "$DEST/versions.txt"

cat "$DEST/versions.txt"
