#!/bin/bash
# T022-A KomegaLegacy solver-entry proof: the Rust cdylib replaces
# KomegaLegacy.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `komegalegacy_` + `komegalegacy_init_`);
#  2. three distinct .sif differentials are byte-identical to the
#     reference except wall-clock lines/values — config A (steady,
#     default legacy bubbles, 1 nonlinear iteration); config B (steady,
#     `Element = "p:1 b:1"` p-bubbles, `Noslip wall BC` + `Omega Wall BC`
#     wall treatment, 2 nonlinear iterations); config C (transient
#     BDF(1) x2 steps, where KOmegaLegacy_Init injects the two-solve
#     minimum); all three carry a SaveScalars solver so the solved
#     Kinetic energy / Kinetic Dissipation fields' min/max land in
#     scalars.dat;
#  3. second full process runs of configs A and C against the port are
#     byte-identical (idempotent re-entry through the SAVE state,
#     including the bx/bxprev bubble history and the BubbleTimestep
#     guard);
#  4. the A/B/C branch variants produce different outputs on both arms
#     (branch-liveness control).
#
# Defined path: all configs run Cartesian coordinates. KomegaLegacy's
# Metric/Symb/dSymb reads sit strictly inside non-Cartesian ELSE
# branches (CoordinateSystemInfo init, SecondInvariant, Metric-weighted
# diffusion), so the Cartesian path reads no uninitialized memory —
# verified by inspection (contrast FDiffusion3D). The port transcribes
# the same branch structure, so the defined path is byte-identical.
#
# Norm class (merged, wall-clock only): wall-clock report lines
# (`Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT`) are
# dropped, and only the numeric timing values after `: (s)` tags
# (FlowSolve per-iteration `Assembly:/Solve: (s)` values) are masked —
# line structure and iteration counts are still compared.
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
WORK="$LANE/scratch/proof-komegalegacy"
CTARGET="${CARGO_TARGET_DIR:-$LANE/scratch/cargo-target-komegalegacy}"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A KomegaLegacy proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: proof pinned to gfortran 13.3.0, re-measure on mismatch"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly (feaadeeac 2026-09-19)"*) ;;
  *) echo "FAIL: unexpected rustc $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" CARGO_TARGET_DIR="$CTARGET" cargo +nightly build --release -q)
PORT_SO="$CTARGET/release/libkomegalegacy.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/KomegaLegacy.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the KomegaLegacy module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T komegalegacy_" "$WORK/port.syms" \
  || { echo "FAIL: komegalegacy_ not the exported entry"; exit 1; }
grep -qx "T komegalegacy_init_" "$WORK/port.syms" \
  || { echo "FAIL: komegalegacy_init_ not exported"; exit 1; }
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
  Steady State Max Iterations = Integer 1
  Output Intervals = Integer 1
  Max Output Level = Integer 4
End
Body 1
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 2
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 3
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Body 4
  Equation = 1
  Material = 1
  Initial Condition = 1
End
Initial Condition 1
  Velocity 1 = Real 0.0
  Velocity 2 = Real 0.0
  Kinetic energy = Real 0.003
  Kinetic Dissipation = Real 0.4
End
Equation 1
  Active Solvers(3) = Integer 1 2 3
End
Material 1
  Density = Real 1.0
  Viscosity = Real 0.01
End
Boundary Condition 1
  Target Boundaries(3) = Integer 1 2 4
  Velocity 1 = Real 0.0
  Velocity 2 = Real 0.0
End
Boundary Condition 2
  Target Boundaries(1) = Integer 3
  Velocity 1 = Real 1.0
  Velocity 2 = Real 0.0
End
Solver 1
  Equation = Navier-Stokes
  Procedure = "FlowSolve" "FlowSolver"
  Stabilization Method = String Stabilized
  Nonlinear System Max Iterations = Integer 20
  Nonlinear System Convergence Tolerance = Real 1.0e-8
  Nonlinear System Newton After Iterations = Integer 3
  Nonlinear System Newton After Tolerance = Real 1.0e-3
  Nonlinear System Relaxation Factor = Real 1.0
  Linear System Solver = Direct
  Linear System Direct Method = UMFPack
  Steady State Convergence Tolerance = Real 1.0e-5
End
Solver 2
  Equation = Komega Legacy
  Procedure = "KomegaLegacy" "KOmegaLegacy"
  Variable = "KE[Kinetic energy:1 Kinetic Dissipation:1]"
  Nonlinear System Max Iterations = Integer 1
  Nonlinear System Convergence Tolerance = Real 1.0e-5
  Linear System Solver = Direct
  Linear System Direct Method = UMFPack
  Steady State Convergence Tolerance = Real 1.0e-4
End
Solver 3
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Velocity 1
  Operator 1 = min
  Variable 2 = Velocity 1
  Operator 2 = max
  Variable 3 = Kinetic energy
  Operator 3 = min
  Variable 4 = Kinetic energy
  Operator 4 = max
  Variable 5 = Kinetic Dissipation
  Operator 5 = min
  Variable 6 = Kinetic Dissipation
  Operator 6 = max
End
EOF
sed -e 's/Nonlinear System Max Iterations = Integer 1/Nonlinear System Max Iterations = Integer 2/' \
    -e 's|Variable = "KE\[Kinetic energy:1 Kinetic Dissipation:1\]"|Variable = "KE[Kinetic energy:1 Kinetic Dissipation:1]"\n  Element = "p:1 b:1"|' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
python3 - "$WORK/caseB.sif" <<'PYEOF'
import sys
p = sys.argv[1]
s = open(p).read()
old = """Boundary Condition 1
  Target Boundaries(3) = Integer 1 2 4
  Velocity 1 = Real 0.0
  Velocity 2 = Real 0.0
End"""
new = """Boundary Condition 1
  Target Boundaries(3) = Integer 1 2 4
  Velocity 1 = Real 0.0
  Velocity 2 = Real 0.0
  Noslip wall BC = Logical True
  Omega Wall BC = Logical True
End"""
assert old in s
s = s.replace(old, new)
open(p, 'w').write(s)
PYEOF
grep -c "Omega Wall BC" "$WORK/caseB.sif"
grep -c 'Element = "p:1 b:1"' "$WORK/caseB.sif"
sed -e 's/Simulation Type = Steady State/Simulation Type = Transient/' \
    -e 's/Steady State Max Iterations = Integer 1/Timestepping Method = BDF\n  BDF Order = Integer 1\n  Timestep Sizes = Real 0.05\n  Timestep Intervals = Integer 2/' \
    -e '/Nonlinear System Max Iterations = Integer 1/d' \
    -e 's/Operator 6 = max/Operator 6 = max\n  Exec Solver = After Timestep/' \
    "$WORK/caseA.sif" > "$WORK/caseC.sif"
grep -c "Timestep Intervals" "$WORK/caseC.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/KomegaLegacy.so"
cp "$PORT_SO" "$WORK/shadow-port/KomegaLegacy.so"
nm -D --defined-only "$WORK/shadow-port/KomegaLegacy.so" | grep -q "T komegalegacy_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }
nm -D --defined-only "$WORK/shadow-port/KomegaLegacy.so" | grep -q "T komegalegacy_init_" \
  || { echo "FAIL: staged port .so lost its init entry"; exit 1; }

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
norm() { # drop wall-clock lines, mask only `: (s)` timing values
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT" "$1" | \
    sed -E 's/(: +\(s\))( +[0-9]+\.[0-9]+)+/\1 TIMING/' || true
}

echo "-- differential A (default legacy bubbles) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-orig/$f" > "$WORK/a0.txt"
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
done
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (p-bubbles + wall BCs, 2 iterations)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- differential C (transient, Init-injected 2 solves) + rerun"
run_case caseC orig runC-orig
run_case caseC port runC-port
run_case caseC port runC-port2
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runC-orig/$f" > "$WORK/c0.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c1.txt"
  diff "$WORK/c0.txt" "$WORK/c1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c0.txt" "$WORK/c1.txt" | head; exit 1; }
done
echo "   config C byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runC-port/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port2/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: transient rerun $f diverges"; exit 1; }
done
echo "   transient rerun byte-identical (bx/BubbleTimestep re-entry)"

echo "-- branch-liveness control (Result Norm differs across A/B/C, both arms)"
for arm in orig port; do
  for cfg in A B C; do
    grep -h "Result Norm   :" "$WORK"/run${cfg}-${arm}/stdout.txt > "$WORK/norm-${cfg}-${arm}.txt"
  done
  diff -q "$WORK/norm-A-${arm}.txt" "$WORK/norm-B-${arm}.txt" > /dev/null 2>&1 && \
    { echo "FAIL: liveness control A-vs-B identical on $arm"; exit 1; }
  diff -q "$WORK/norm-A-${arm}.txt" "$WORK/norm-C-${arm}.txt" > /dev/null 2>&1 && \
    { echo "FAIL: liveness control A-vs-C identical on $arm"; exit 1; }
done
echo "   control holds: branches change the numerics in both arms"

echo "KOMEGALEGACY PROOF GREEN wall=$((SECONDS - T0))s"
