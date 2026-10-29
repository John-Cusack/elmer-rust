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
for tool in cc gcc as ar ld ln cp rm mkdir mv head grep ldd nm find ls dirname basename cargo rustc python3 mpicc; do
    p="$(command -v "$tool" || true)"
    if [ -n "$p" ]; then ln -s "$p" "$NOFCDIR/$tool"; fi
done
export PATH="$NOFCDIR"
# Hermetic build: drop ambient compiler wrappers/caches.
unset RUSTC_WRAPPER SCCACHE_DIR SCCACHE_CACHE_SIZE CARGO_INCREMENTAL
# Only workspace-member crates (lanes merge independently; a hardcoded
# list breaks every lane whose crates have not merged yet, and excluded
# crates such as elmer-oracle are not buildable with -p).
PKGS=""
for name in $(cargo metadata --no-deps --format-version 1 \
    --manifest-path "$MANIFEST" 2>/dev/null | \
    python3 -c 'import json,sys; print(" ".join(p["name"] for p in json.load(sys.stdin)["packages"]))'); do
    PKGS="$PKGS -p $name"
done
[ -n "$PKGS" ] || fail "no native workspace members found"
# shellcheck disable=SC2086
CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}" \
    cargo build --manifest-path "$MANIFEST" $PKGS \
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
