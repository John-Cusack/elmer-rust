#!/bin/sh
# Zero-Fortran enforcement for the Rust-native workspace (S1 front-end owns it).
# Fails when any native crate needs gfortran/libgfortran to build or link.
# Usage: sh native/check-zero-fortran.sh   (run from the worktree root)
set -eu

ROOT="$(dirname "$0")"
MANIFEST="$ROOT/Cargo.toml"
TARGET="$ROOT/target"

fail() { echo "zero-fortran FAIL: $1" >&2; exit 1; }

echo "== 1. dependency scan: no crate may pull a Fortran toolchain =="
# `cc`/`cmake` build scripts are the usual Fortran-smuggling path; the
# front-end vendors C sources and compiles them with a direct `cc` call.
if grep -Ei '^\s*cc\s*=|^\s*cmake\s*=' "$ROOT"/elmer-*/Cargo.toml 2>/dev/null; then
    fail "build-dependency on cc/cmake found (may invoke a Fortran compiler)"
fi
echo "ok: no cc/cmake build dependencies"

echo "== 2. build with gfortran absent from PATH =="
NOFCDIR="$(mktemp -d)"
trap 'rm -rf "$NOFCDIR"' EXIT INT TERM
for tool in cc gcc ar ld ln cp rm mkdir mv head grep ldd nm find ls cargo rustc; do
    p="$(command -v "$tool" || true)"
    if [ -n "$p" ]; then ln -s "$p" "$NOFCDIR/$tool"; fi
done
# PATH contains ONLY the stub dir: gfortran must NOT resolve.
export PATH="$NOFCDIR"
# Belt and braces: any $(FC) use explodes loudly.
export FC=false
command -v gfortran >/dev/null 2>&1 && fail "gfortran still on PATH during check"
CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}" \
    cargo build --manifest-path "$MANIFEST" -p elmer-sif -p elmer-mesh -p elmer-fe \
    || fail "native build needs something outside the gfortran-free PATH"
echo "== 3. linked-binary scan: no libgfortran =="
EXAMPLE_BIN="$(ls "$TARGET"/debug/examples/sifproof 2>/dev/null || true)"
if [ -n "$EXAMPLE_BIN" ]; then
    if ldd "$EXAMPLE_BIN" 2>/dev/null | grep -qi fortran; then
        ldd "$EXAMPLE_BIN" | grep -i fortran >&2
        fail "libgfortran linked into $EXAMPLE_BIN"
    fi
    echo "ok: $EXAMPLE_BIN links no Fortran runtime"
else
    echo "skip: example binary not present (build it to enable the ldd scan)"
fi
echo "== 4. staticlib symbol scan =="
MATC_A="$(find "$TARGET"/debug/build -name libmatc.a 2>/dev/null | head -1 || true)"
if [ -n "$MATC_A" ]; then
    if nm "$MATC_A" 2>/dev/null | grep -qiE 'gfortran|_gfortran'; then
        fail "gfortran symbols in $MATC_A"
    fi
    echo "ok: no gfortran symbols in vendored libmatc.a"
else
    echo "skip: libmatc.a not present"
fi
echo "zero-fortran PASS"
