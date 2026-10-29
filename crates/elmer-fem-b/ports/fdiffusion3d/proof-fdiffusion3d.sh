#!/bin/bash
# T022-A FDiffusion3D solver-entry proof: the Rust cdylib replaces FDiffusion3D.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `fourierdiffusion3dsolver_`);
#  2. two distinct clean-room .sif differentials on a REAL 3-D mesh
#     (config A: Angular Frequency + Electric Conductivity + 6-component
#     body force; config B: Frequency branch + Electrical Conductivity
#     Warn branch + new values + 2 nonlinear iterations) are
#     byte-identical to the reference except wall-clock lines and
#     per-iteration (s) timing values;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. the B config produces different output than A on both arms
#     (branch-liveness control).
#
# REFERENCE-UB DOCUMENTATION (per the Wave-1 "defined path" rule):
# FDiffusion3D.F90's LocalMatrix initializes Metric to identity but never
# calls CoordinateSystemInfo on the Cartesian path, so in Cartesian
# coordinates it reads uninitialized stack arrays Symb/dSymb in every
# element assembly. The reference is therefore nondeterministic across
# processes in Cartesian (and no port can match it there — the same binary
# disagrees with itself run to run). The author's evident intent is
# unambiguous: CoordinateSystemInfo's Cartesian branch sets Symbols = 0,
# dSymbols = 0, so the wrapper zeroes both next to the Metric init
# (bit-identical with the reference wherever the reference IS defined:
# every non-Cartesian path calls CoordinateSystemInfo, which overwrites
# all four). Accordingly both differential configs run the defined path
# `Coordinate System = Cylindrical` (the 3-D cylindric system) on the
# unit cube shifted to x in [1,2] (r never 0, so the guarded 1/r terms
# never engage); that path exercises the full Metric/Symb tensor coupling.
#
# The raw `Assembly: N % done` progress WRITE is gated on RealTime()-at0 >
# 1.0 s, so its presence and percentage are pure wall-clock artifacts
# (the 2-D sibling's equivalent goes through Info level 5 and is
# suppressed by Max Output Level = 4; here the original writes it raw, so
# norm() drops those lines — wall-clock class only, nothing else).
#
# The campaign Procedure form (`File "FDiffusion3D"
# "FourierDiffusion3DSolver"`, library-vs-entry names differing) stages the
# port .so under its own reference name FDiffusion3D.so.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent)
# and the Elmer checkout only for the cube.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-fdiffusion3d"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/cube*
echo "== T022-A FDiffusion3D proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: need gfortran 13.3.0, got: $GFORTRAN_VERSION"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: need rustc 1.100.0-nightly, got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$PORT_DIR/target}" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libfdiffusion3d.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/FDiffusion3D.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the FDiffusion3D module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T fourierdiffusion3dsolver_" "$WORK/port.syms" \
  || { echo "FAIL: fourierdiffusion3dsolver_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's cube.grd, shifted to x in [1,2])"
cp "$ELMER_FORK/fem/tests/radiator3d/cube.grd" "$WORK/"
sed -i 's/^Subcell Limits 1 = 0 1$/Subcell Limits 1 = 1 2/' "$WORK/cube.grd"
grep -q "^Subcell Limits 1 = 1 2$" "$WORK/cube.grd" || { echo "FAIL: mesh shift did not apply"; exit 1; }
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 cube > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "cube"
End
Simulation
  Coordinate System = Cylindrical
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
  Max Output Level = Integer 4
End
Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(2) = Integer 1 2
  FourierDiffusion3D = Logical True
  Angular Frequency = Real 100.0
End
Material 1
  Electric Conductivity = Real 5.0e6
  Magnetic Permeability = Real 1.2566e-6
End
Body Force 1
  AC Magnetic Field 1 Re = Real 1.0
  AC Magnetic Field 1 Im = Real 0.5
  AC Magnetic Field 2 Re = Real 0.25
  AC Magnetic Field 2 Im = Real -0.5
  AC Magnetic Field 3 Re = Real 0.75
  AC Magnetic Field 3 Im = Real 0.125
End
Solver 1
  Equation = FourierDiffusion3D
  Procedure = File "FDiffusion3D" "FourierDiffusion3DSolver"
  Variable = BField
  Variable DOFs = 6
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = Real 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = BField
  Operator 1 = min
  Variable 2 = BField
  Operator 2 = max
End
EOF
# Config B: Frequency (not Angular Frequency) keyword branch, Electrical
# Conductivity (Warn) branch, new body-force values, and 2 nonlinear
# iterations (SAVE re-entry inside one process).
sed -e 's/Angular Frequency = Real 100.0/Frequency = Real 15.0/' \
    -e 's/Electric Conductivity = Real 5.0e6/Electrical Conductivity = Real 3.0e6/' \
    -e 's/Magnetic Permeability = Real 1.2566e-6/Magnetic Permeability = Real 2.0e-6/' \
    -e 's/AC Magnetic Field 1 Re = Real 1.0/AC Magnetic Field 1 Re = Real 2.0/' \
    -e 's/AC Magnetic Field 1 Im = Real 0.5/AC Magnetic Field 1 Im = Real -1.0/' \
    -e 's/AC Magnetic Field 2 Re = Real 0.25/AC Magnetic Field 2 Re = Real 0.5/' \
    -e 's/AC Magnetic Field 2 Im = Real -0.5/AC Magnetic Field 2 Im = Real 1.0/' \
    -e 's/AC Magnetic Field 3 Re = Real 0.75/AC Magnetic Field 3 Re = Real 1.5/' \
    -e 's/AC Magnetic Field 3 Im = Real 0.125/AC Magnetic Field 3 Im = Real -0.25/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
python3 - "$WORK/caseB.sif" <<'PYEOF'
import sys
p = sys.argv[1]
t = open(p).read()
anchor = "  Linear System Preconditioning = ILU0\nEnd"
assert anchor in t
t = t.replace(anchor, "  Linear System Preconditioning = ILU0\n  Nonlinear System Max Iterations = Integer 2\nEnd", 1)
open(p, "w").write(t)
PYEOF
grep -c "Nonlinear System Max Iterations" "$WORK/caseB.sif"
# Config B keeps both conductivity spellings: 'Electrical' fires the
# Warn branch, then (per the 3-D unit) 'Electric' overwrites the value.
python3 - "$WORK/caseB.sif" <<'PYEOF'
import sys
p = sys.argv[1]
t = open(p).read()
anchor = "  Electrical Conductivity = Real 3.0e6\n"
assert anchor in t
t = t.replace(anchor, anchor + "  Electric Conductivity = Real 4.0e6\n", 1)
open(p, "w").write(t)
PYEOF
grep -c "Electric Conductivity" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/FDiffusion3D.so"
cp "$PORT_SO" "$WORK/shadow-port/FDiffusion3D.so"
nm -D --defined-only "$WORK/shadow-port/FDiffusion3D.so" | grep -q "T fourierdiffusion3dsolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/cube" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines, mask per-iteration (s) timing values
  # Wall-clock class only: STARTED/FINISHED/TOTAL lines, the (s) values on
  # the deterministic iter/timing structure lines, and the raw RealTime-gated
  # `Assembly: N % done` progress lines (presence and percentage are pure
  # timing artifacts of the >1 s gate).
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" \
    | grep -v -E "^   Assembly: +[0-9]+ % done$" \
    | sed -E '/iter:.*(Assembly|Solve): +\(s\)/ s/[0-9]+\.[0-9]+/T/g' || true
}

echo "-- differential A (Angular Frequency branch) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stderr.txt scalars.dat; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 \
    || { echo "FAIL: config A $f diverges"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 \
  || { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
norm "$WORK/runA-port2/stdout.txt" > "$WORK/a2.txt"
diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 \
  || { echo "FAIL: rerun stdout diverges"; diff "$WORK/a1.txt" "$WORK/a2.txt" | head; exit 1; }
for f in stderr.txt scalars.dat; do
  diff "$WORK/runA-port/$f" "$WORK/runA-port2/$f" > /dev/null 2>&1 \
    || { echo "FAIL: rerun $f diverges"; diff "$WORK/runA-port/$f" "$WORK/runA-port2/$f" | head; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (Frequency + Warn + 2-iteration branches)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  if [ "$f" = stdout.txt ]; then
    norm "$WORK/runB-orig/$f" > "$WORK/b0.txt"
    norm "$WORK/runB-port/$f" > "$WORK/b1.txt"
    diff "$WORK/b0.txt" "$WORK/b1.txt" > /dev/null 2>&1 \
      || { echo "FAIL: config B $f diverges"; diff "$WORK/b0.txt" "$WORK/b1.txt" | head; exit 1; }
  else
    diff "$WORK/runB-orig/$f" "$WORK/runB-port/$f" > /dev/null 2>&1 \
      || { echo "FAIL: config B $f diverges"; diff "$WORK/runB-orig/$f" "$WORK/runB-port/$f" | head; exit 1; }
  fi
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- liveness control (B output differs from A on both arms)"
if diff -q "$WORK/runA-orig/scalars.dat" "$WORK/runB-orig/scalars.dat" > /dev/null; then
  echo "FAIL: liveness control dead (orig A == orig B)"; exit 1;
fi
if diff -q "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" > /dev/null; then
  echo "FAIL: liveness control dead (port A == port B)"; exit 1;
fi
echo "   liveness control live on both arms"

T1=$SECONDS
echo "FDIFFUSION3D PROOF GREEN (wall: $((T1-T0))s)"
