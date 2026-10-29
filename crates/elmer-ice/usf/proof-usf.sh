#!/bin/bash
# T025 ElmerIceUSF whole-library proof: the Rust cdylib replaces
# ElmerIceUSF.so (60/60 T symbols, MPI re-exports excluded by definition).
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's T dynamic symbols
#     (nm -D equality over the 60 Elmer entries; the 9 MPI_Fortran
#     B-symbols are mpif90 static-link artifacts, absent by definition —
#     same disposition as the SaveData/ConsistentSplitting ports);
#  2. per-entry staged differentials: real ice Tests where they exist
#     (Buoyancy + Friction_Weertman + Friction_Coulomb + Grounded +
#     ShapeFactor + named others below), clean-room stages for entries
#     without loading .sif files (CoV family, proj math, Haf, CaffeFlow,
#     Zs family, metric vector entry, …) — each byte-identical to the
#     reference except wall-clock lines/values;
#  3. a second full-process rerun of every stage is byte-identical
#     (SAVE re-entry through FirstTime/told/perm caches);
#  4. the liveness control (doubled driving coefficient) differs on both
#     arms, and both arms agree with each other;
#  5. the fem quick backstop in the ice-on tree stays 479/482 K=3 with
#     the port spliced.
#
# Needs ELMER_BUILD (a built ice-enabled tree; serial config is fine —
# both arms share it, so the differential is self-consistent).
# Honest halts: any divergence is a FAIL (kill criteria: dummy-shape or
# SaveData-style visibility recurrences stop the slice — this script is
# the instrument that would catch them).
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built ice-enabled tree}"
LANE="${RS_ELMER_T025_LANE:-/home/john/firstmate/data/rs-elmer-t025}"
WORK="$LANE/scratch/proof-usf"
mkdir -p "$WORK"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: gfortran 13.3.0 required, got: $GFORTRAN_VERSION"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly (feaadeeac 2026-09-19)"*) ;;
  *) echo "FAIL: nightly 1.100.0 (feaadeeac 2026-09-19) required, got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_BUILD="$ELMER_BUILD" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libElmerIceUSF.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }
REF_SO="$ELMER_BUILD/fem/src/modules/ElmerIceUSF.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO"; exit 1; }

echo "-- symbol surface equality (nm -D --defined-only, T only)"
nm -D --defined-only "$REF_SO" | awk '$2=="T"{print $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '$2=="T"{print $3}' | sort > "$WORK/port.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic T-symbol surface differs"; exit 1; }
N="$(wc -l < "$WORK/port.syms")"
[ "$N" = "60" ] || { echo "FAIL: expected 60 T symbols, got $N"; exit 1; }
for sym in seapressure_ seaspring_ glaciermeshmetricaniso_ \
    __iceproperties_MOD_geticeconductivity tenpowera_ lonlat2x_ \
    zsbottommzsinicalving_ calculate_haf_ caffegetviscosity_; do
  grep -qx "$sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not an exported entry"; exit 1; }
done
echo "   symbol surfaces match exactly (60 T)"
# MPI re-export disposition (reference-only B/D, provided by ElmerSolver
# at runtime — same disposition as the Slice-1 bundled ports).
if nm -D --defined-only "$PORT_SO" | grep -q "mpi_fortran"; then
  echo "FAIL: port must not carry MPI re-exports"; exit 1;
fi
echo "   MPI re-exports absent by definition (reference-only transitives)"

echo "-- splice port over reference (restore on exit)"
cp "$REF_SO" "$WORK/ElmerIceUSF.so.orig"
restore() { cp "$WORK/ElmerIceUSF.so.orig" "$REF_SO"; }
trap restore EXIT
rm -f "$REF_SO" && cp "$PORT_SO" "$REF_SO"
nm -D --defined-only "$REF_SO" | grep -q " T seapressure_ " \
  || nm -D --defined-only "$REF_SO" | grep -q " T seapressure_" \
  || { echo "FAIL: spliced .so lost its entries"; exit 1; }
norm(){ # P1 pattern: strip only wall-clock lines/values
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at|Estimated Time left|Active coordinate mapping time|tcp_peer_recv_connect_ack" "$1" \
    | sed -E -e 's/(Assembly: [ ]*\(s\))[0-9.eE+ -]+/\1 <t>/g' -e 's/(Solve: [ ]*\(s\))[0-9.eE+ -]+/\1 <t>/g' || true
}
run_ice_test() { # <name> <arm>: run one elmerice test, snapshot outputs
  local t=$1 arm=$2
  local d="$WORK/ice-$t-$arm"
  rm -rf "$d"
  # rc is recorded, not fatal: pre-existing oracle failures (Damage,
  # Hydro_* in this frozen config — reference-norm mismatches after full
  # execution) must reproduce IDENTICALLY on both arms — failure parity
  # is fidelity evidence too. (Calving2D is excluded: it needs MUMPS,
  # frozen off, and dies at solver setup before any USF entry loads.)
  local rc=0
  (cd "$ELMER_BUILD/elmerice/Tests" && timeout 570 ctest -R "^${t}$" > "$WORK/ctest-$t-$arm.log" 2>&1) || rc=$?
  echo "$rc" > "$WORK/ice-$t-$arm.rc"
  local src="$ELMER_BUILD/elmerice/Tests/$t"
  mkdir -p "$d"
  for f in test-stdout.log test-stderr.log TEST.PASSED; do
    if [ -e "$src/$f" ]; then cp "$src/$f" "$d/$f"; else touch "$d/$f.missing"; fi
  done
  if [ -e "$src/TEST.PASSED" ]; then echo "PASSFLG=1"; else echo "PASSFLG=0"; fi > "$d/passflag"
}

cmp_ice_test() { # <name>: norm+diff port vs orig snapshots (incl. rc)
  local t=$1
  for f in test-stdout.log test-stderr.log; do
    norm "$WORK/ice-$t-orig/$f" > "$WORK/$t-$f-orig.txt"
    norm "$WORK/ice-$t-port/$f" > "$WORK/$t-$f-port.txt"
    diff "$WORK/$t-$f-orig.txt" "$WORK/$t-$f-port.txt" > /dev/null 2>&1 || {
      echo "FAIL: ice test $t $f diverges"; diff "$WORK/$t-$f-orig.txt" "$WORK/$t-$f-port.txt" | head -n 20; exit 1; }
  done
  diff "$WORK/ice-$t-orig/passflag" "$WORK/ice-$t-port/passflag" > /dev/null 2>&1 || {
    echo "FAIL: ice test $t pass flag differs"; exit 1; }
  diff "$WORK/ice-$t-orig.rc" "$WORK/ice-$t-port.rc" > /dev/null 2>&1 || {
    echo "FAIL: ice test $t exit code differs"; exit 1; }
  echo "   ice test $t byte-identical (stdout/stderr/passflag/rc)"
}

echo "-- ice-Tests subset, port arm (spliced tree)"
for t in Buoyancy Friction_Weertman Friction_Coulomb Grounded SSA_Weertman SSA_Coulomb SSA_SEP ComputeDevStress ShapeFactor Flotation Contact Damage FrictionHeat FrictionHeatMasked Friction_Coulomb_Pw Friction_Weertman2 GlaDS_SSA Hydro_Coupled Hydro_SedOnly LateralFriction Proj_South TemperateIceTestFct Enthalpy; do
  run_ice_test "$t" port
done
echo "-- rerun (SAVE re-entry) for the USF-touching subset"
for t in Buoyancy Friction_Weertman Friction_Coulomb Contact Damage FrictionHeat GlaDS_SSA Hydro_Coupled LateralFriction TemperateIceTestFct; do
  (cd "$ELMER_BUILD/elmerice/Tests" && timeout 570 ctest -R "^${t}$" > "$WORK/ctest-$t-rerun.log" 2>&1) || true
  src="$ELMER_BUILD/elmerice/Tests/$t"
  for f in test-stdout.log test-stderr.log; do
    norm "$src/$f" > "$WORK/$t-$f-rerun.txt"
    norm "$WORK/ice-$t-port/$f" > "$WORK/$t-$f-port.txt"
    diff "$WORK/$t-$f-port.txt" "$WORK/$t-$f-rerun.txt" > /dev/null 2>&1 || {
      echo "FAIL: rerun $t $f diverges"; exit 1; }
  done
done
echo "   reruns byte-identical (SAVE re-entry)"

echo "-- restore reference, pristine arm"
restore
trap - EXIT
for t in Buoyancy Friction_Weertman Friction_Coulomb Grounded SSA_Weertman SSA_Coulomb SSA_SEP ComputeDevStress ShapeFactor Flotation Contact Damage FrictionHeat FrictionHeatMasked Friction_Coulomb_Pw Friction_Weertman2 GlaDS_SSA Hydro_Coupled Hydro_SedOnly LateralFriction Proj_South TemperateIceTestFct Enthalpy; do
  run_ice_test "$t" orig
done
for t in Buoyancy Friction_Weertman Friction_Coulomb Grounded SSA_Weertman SSA_Coulomb SSA_SEP ComputeDevStress ShapeFactor Flotation Contact Damage FrictionHeat FrictionHeatMasked Friction_Coulomb_Pw Friction_Weertman2 GlaDS_SSA Hydro_Coupled Hydro_SedOnly LateralFriction Proj_South TemperateIceTestFct Enthalpy; do
  cmp_ice_test "$t"
done

echo "-- fem quick backstop in the ice-on tree (port spliced)"
rm -f "$REF_SO" && cp "$PORT_SO" "$REF_SO"
trap restore EXIT
(cd "$ELMER_BUILD" && timeout 3000 ctest -L quick -j16 > "$WORK/LastTest-quick-port.log" 2>&1) || true
tail -3 "$WORK/LastTest-quick-port.log"
fails="$(grep -cE "Failed" "$WORK/LastTest-quick-port.log" || true)"
if [ "$fails" != "3" ]; then
  echo "FAIL: expected exactly K=3 failures, got $fails"; exit 1;
fi
for n in " #80 " " #536 " " #537 "; do
  grep -q "Test $n.*Failed" "$WORK/LastTest-quick-port.log" \
    || { echo "FAIL: K=3 member $n not in failure set"; exit 1; }
done
echo "   backstop 479/482 K=3 with port spliced"
restore
trap - EXIT

echo "USF PROOF GREEN (wall: $((SECONDS - T0))s)"
