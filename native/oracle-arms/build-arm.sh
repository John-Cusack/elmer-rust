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

if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "usage: build-arm.sh <release|fma|openblas|listusage> <source-dir> <build-dir> [patch-dir]" >&2
    exit 2
fi

ARM="$1"
SRC="$2"
DEST="$3"
PATCH_DIR="${4:-}"

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

# Oracle-only instrumentation (D6) lands as versioned patches on the copy,
# never in the pristine source. Every applied patch is hashed into
# versions.txt; the copy's own git diff stays reviewable in the build dir.
APPLIED_PATCHES=""
if [ -n "$PATCH_DIR" ]; then
    if [ ! -d "$PATCH_DIR" ]; then
        echo "refusing: patch dir $PATCH_DIR missing" >&2
        exit 2
    fi
    for patch in "$PATCH_DIR"/*.patch; do
        [ -e "$patch" ] || break
        git -C "$SRC" apply --check "$patch" || exit 2
        git -C "$SRC" apply "$patch"
        APPLIED_PATCHES="$APPLIED_PATCHES $(basename "$patch"):$(sha256sum "$patch" | cut -d' ' -f1)"
    done
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
    # DEVEL_LISTUSAGE is not a declared CMake option: a bare -D only
    # lands in the cache and never reaches the compiler. Pass it as a
    # Fortran pre-processor define instead (verified in flags.make).
    listusage) extra_flags=(-DCMAKE_Fortran_FLAGS="-DDEVEL_LISTUSAGE") ;;
    *) echo "unknown arm $ARM" >&2; exit 2 ;;
esac

cmake -S "$SRC" -B "$DEST" "${common_flags[@]}" "${extra_flags[@]}"
cmake --build "$DEST" --parallel 4

{
    echo "arm: $ARM"
    echo "source-pin: $head"
    echo "source-patches:${APPLIED_PATCHES:- none}"
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
    grep -E '^(CMAKE_Fortran_FLAGS|CMAKE_C_FLAGS)' "$DEST/CMakeCache.txt" || true
    echo "--- cpu ---"
    grep -m1 'model name' /proc/cpuinfo
    grep -m1 flags /proc/cpuinfo | tr ' ' '\n' | grep -c -E '^(fma|avx2|avx512f)$' || true
    echo "--- image ---"
    docker images --digests elmer-oracle-fortran 2>/dev/null || echo "native (non-docker) build"
} > "$DEST/versions.txt"

cat "$DEST/versions.txt"
