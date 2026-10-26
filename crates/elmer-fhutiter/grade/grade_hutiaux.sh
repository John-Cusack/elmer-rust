#!/bin/sh
# T010 driver-level grade for `hutiaux` (lane-fhutiter).
# Substitutes the port for the reference `huti_aux.F90.o` at link time
# (the tool's `substitute` halt on Fortran CHECK-format DependInfo is a
# recorded tool gap — mechanics here are the identical `ld -r` composition),
# then diffs stdout+stderr+exit over the deterministic corpus.
#
# Layout (all paths relative to the elmer-rust checkout root unless absolute):
#   $SRC   pristine Elmer tree @ 9f6af2f85 (read-only)
#   $BUILD configured cmake build of $SRC (provides huti_aux.F90.o + .mod)
#   $CRATE this lane's ports/hutiaux
# Outputs: grade/probe_ref.{out,err}, grade/probe_port.{out,err}, verdict.
set -eu
: "${SRC:?set SRC to the pristine Elmer tree}"
: "${BUILD:?set BUILD to the configured cmake build of SRC}"
LANE="$(cd "$(dirname "$0")/.." && pwd)"
CRATE="$LANE/ports/hutiaux"
GRADE="$LANE/grade"
MODDIR="$BUILD/fmodules"
REFOBJ="$BUILD/fhutiter/src/CMakeFiles/fhuti.dir/huti_aux.F90.o"
LIB="$LANE/target/release/libhutiaux.a"

"$LANE/build-check.sh"
test -f "$REFOBJ" || { echo "missing reference object: $REFOBJ" >&2; exit 1; }

echo "--- build port (pinned nightly) ---"
/usr/bin/time -v cargo build --offline --release --manifest-path "$LANE/Cargo.toml" 2> "$GRADE/port-build.time"
test -f "$LIB" || { echo "missing port archive: $LIB" >&2; exit 1; }

echo "--- compile probe ---"
gfortran -O0 -I "$MODDIR" -c "$GRADE/probe_hutiaux.F90" -o "$GRADE/probe.o"

echo "--- link reference (orig object) ---"
gfortran -O0 -o "$GRADE/probe_ref" "$GRADE/probe.o" "$REFOBJ" -lgfortran

echo "--- link port (whole-archive staticlib) ---"
gfortran -O0 -o "$GRADE/probe_port" "$GRADE/probe.o" \
  -Wl,--whole-archive "$LIB" -Wl,--no-whole-archive -lgfortran

echo "--- run both (unbuffered: true program order across languages) ---"
export GFORTRAN_UNBUFFERED_ALL=Y
/usr/bin/time -v -o "$GRADE/probe_ref.time" "$GRADE/probe_ref" > "$GRADE/probe_ref.out" 2> "$GRADE/probe_ref.err"; echo $? > "$GRADE/probe_ref.rc"
/usr/bin/time -v -o "$GRADE/probe_port.time" "$GRADE/probe_port" > "$GRADE/probe_port.out" 2> "$GRADE/probe_port.err"; echo $? > "$GRADE/probe_port.rc"

echo "--- verdict ---"
fail=0
cmp "$GRADE/probe_ref.out" "$GRADE/probe_port.out" || { echo "STDOUT DIFFERS"; fail=1; }
cmp "$GRADE/probe_ref.err" "$GRADE/probe_port.err" || { echo "STDERR DIFFERS"; fail=1; }
cmp "$GRADE/probe_ref.rc" "$GRADE/probe_port.rc" || { echo "EXIT DIFFERS"; fail=1; }

echo "--- strong-symbol gate ---"
# (a) splice completeness: port defines exactly the 12 removed unit symbols.
nm -g --defined-only "$REFOBJ" | awk '$2=="T"{print $3}' | sort > "$GRADE/ref-T.txt"
nm -g --defined-only "$LIB" | awk '$2=="T" && $3 ~ /^__huti_aux_MOD_/{print $3}' | sort -u > "$GRADE/port-T.txt"
if cmp -s "$GRADE/ref-T.txt" "$GRADE/port-T.txt"; then
  echo "symbol gate: 12/12 unit symbols, zero extras in module namespace"
else
  echo "SYMBOL GATE FAILED"; diff "$GRADE/ref-T.txt" "$GRADE/port-T.txt" || true; fail=1
fi
# (b) user-namespace cleanliness: every other port T symbol must be a
# toolchain class (compiler builtins, rust_eh_personality, alloc shims,
# monomorphized std/gimli internals) resolved at merge by the runbook
# one-provider/byte-identity rule — reported, never asserted away.
nm -g --defined-only "$LIB" | awk '$2=="T"{print $3}' | grep -v -E '^(__huti_aux_MOD_|rust_eh_personality$|__rust_|_R|_GLOBAL_|__absv|__addv|__ashl|__ashr|__bswap|__clz|__cmp|__ctz|__div|__ffsti2|__fix|__float|__lshr|__mod|__mulv|__neg|__parity|__popcount|__subv|__trampoline|__ucmp|__udiv|__umod|__main$|__stack_chk|__gmon_start__|_ITM_)' \
  | sort -u > "$GRADE/port-T-nontoolchain.txt" || true
if [ -s "$GRADE/port-T-nontoolchain.txt" ]; then
  echo "NON-TOOLCHAIN port T symbols:"; cat "$GRADE/port-T-nontoolchain.txt"; fail=1
else
  echo "user namespace clean: all other T symbols are toolchain class"
fi
nm -g --defined-only "$LIB" | awk '$2=="T"' | wc -l | xargs echo "port T total (unit 12 + toolchain):"

if [ "$fail" = 0 ]; then echo "GRADE GREEN"; else echo "GRADE RED"; exit 1; fi
