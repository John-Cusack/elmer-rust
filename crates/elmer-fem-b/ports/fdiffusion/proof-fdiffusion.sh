#!/bin/bash
# T022-A FDiffusion solver-entry proof: the Rust cdylib replaces FDiffusion.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `fourierdiffusionsolver_`);
#  2. two distinct .sif differentials (config A: Angular Frequency +
#     Electric Conductivity + impedance BCs; config B: Frequency branch +
#     Electrical Conductivity Warn branch + flux/impedance values + 2
#     nonlinear iterations) are byte-identical to the reference except
#     wall-clock lines and per-iteration (s) timing values;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. the B config produces different output than A on both arms
#     (branch-liveness control).
#
# The campaign Procedure form (`File "FDiffusion"
# "FourierDiffusionSolver"`, library-vs-entry names differing) stages the
# port .so under its own reference name FDiffusion.so.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent)
# and the Elmer checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-fdiffusion"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A FDiffusion proof: work=$WORK build=$ELMER_BUILD"
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
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: nightly 1.100.0 required, got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$PORT_DIR/target}" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libfdiffusion.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/FDiffusion.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the FDiffusion module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T fourierdiffusionsolver_" "$WORK/port.syms" \
  || { echo "FAIL: fourierdiffusionsolver_ not the exported entry"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's square.grd)"
cp "$ELMER_FORK/fem/tests/PoissonDB/square.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "square"
End
Simulation
  Coordinate System = Cartesian
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
Body 2
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 3
  Equation = 1
  Material = 1
  Body Force = 1
End
Body 4
  Equation = 1
  Material = 1
  Body Force = 1
End
Equation 1
  Active Solvers(2) = Integer 1 2
  FourierDiffusion = Logical True
  Angular Frequency = Real 100.0
End
Material 1
  Electric Conductivity = Real 5.0e6
  Magnetic Permeability = Real 1.2566e-6
End
Body Force 1
  AC Magnetic Field Re = Real 1.0
  AC Magnetic Field Im = Real 0.5
End
Boundary Condition 1
  Target Boundaries(1) = Integer 1
  Wave Impedance 1 = Real 1.0e-3
  Wave Impedance 2 = Real 2.0e-3
  Wave Flux 1 = Real 0.1
  Wave Flux 2 = Real 0.2
End
Boundary Condition 2
  Target Boundaries(3) = Integer 2 3 4
  Wave Impedance 1 = Real 3.0e-3
  Wave Impedance 2 = Real 4.0e-3
  Wave Flux 1 = Real 0.3
  Wave Flux 2 = Real 0.4
End
Solver 1
  Equation = FourierDiffusion
  Procedure = File "FDiffusion" "FourierDiffusionSolver"
  Variable = BField
  Variable DOFs = 2
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
# Conductivity (Warn) branch, new flux/impedance/body-force values, and 2
# nonlinear iterations (SAVE re-entry inside one process).
sed -e 's/Angular Frequency = Real 100.0/Frequency = Real 15.0/' \
    -e 's/Electric Conductivity = Real 5.0e6/Electrical Conductivity = Real 3.0e6/' \
    -e 's/Magnetic Permeability = Real 1.2566e-6/Magnetic Permeability = Real 2.0e-6/' \
    -e 's/AC Magnetic Field Re = Real 1.0/AC Magnetic Field Re = Real 2.0/' \
    -e 's/AC Magnetic Field Im = Real 0.5/AC Magnetic Field Im = Real -1.0/' \
    -e 's/Wave Impedance 1 = Real 1.0e-3/Wave Impedance 1 = Real 5.0e-3/' \
    -e 's/Wave Impedance 2 = Real 2.0e-3/Wave Impedance 2 = Real 6.0e-3/' \
    -e 's/Wave Flux 1 = Real 0.1/Wave Flux 1 = Real 0.7/' \
    -e 's/Wave Flux 2 = Real 0.2/Wave Flux 2 = Real 0.8/' \
    -e 's/Wave Impedance 1 = Real 3.0e-3/Wave Impedance 1 = Real 7.0e-3/' \
    -e 's/Wave Impedance 2 = Real 4.0e-3/Wave Impedance 2 = Real 8.0e-3/' \
    -e 's/Wave Flux 1 = Real 0.3/Wave Flux 1 = Real 0.9/' \
    -e 's/Wave Flux 2 = Real 0.4/Wave Flux 2 = Real 1.0/' \
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

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/FDiffusion.so"
cp "$PORT_SO" "$WORK/shadow-port/FDiffusion.so"
nm -D --defined-only "$WORK/shadow-port/FDiffusion.so" | grep -q "T fourierdiffusionsolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/square" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines, mask per-iteration (s) timing values
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|STARTED AT|File started at" "$1" \
    | sed -E '/iter:.*(Assembly|Solve): \(s\)/ s/[0-9]+\.[0-9]+/T/g' || true
}

echo "-- differential A (Angular Frequency branch) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stderr.txt scalars.dat; do
  diff -q "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 \
    || { echo "FAIL: config A $f differs"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 \
  || { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-port/$f" > "$WORK/r1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/r2.txt"
  diff -q "$WORK/r1.txt" "$WORK/r2.txt" > /dev/null \
    || { echo "FAIL: rerun $f differs (SAVE re-entry)"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (Frequency + Warn + 2-iteration branches)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b0.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b1.txt"
  diff -q "$WORK/b0.txt" "$WORK/b1.txt" > /dev/null \
    || { echo "FAIL: config B $f diverges"; diff "$WORK/b0.txt" "$WORK/b1.txt" | head; exit 1; }
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
echo "FDIFFUSION PROOF GREEN (wall: $((T1-T0))s)"
