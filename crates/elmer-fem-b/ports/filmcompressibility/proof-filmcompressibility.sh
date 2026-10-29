#!/bin/bash
# T022-A FilmCompressibility solver-entry proof: the Rust cdylib replaces
# FilmCompressibility.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `filmcompressibility_` +
#     `filmcompressibility_init_`);
#  2. two distinct .sif differentials are byte-identical to the
#     reference except wall-clock lines — config A (FullMode with AC
#     Diffusion: the PDE-assembly path) and config B (Differential Mode
#     without diffusion plus asymmetric AC relaxation: the diagonal
#     assembly + WHERE-update paths, over 2 BDF timesteps); each carries
#     a `SaveScalars` solver so the solved `AC` field's min/max land in
#     scalars.dat;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state).
# Branch-liveness control: the `AC Average` report and the saved
# scalars differ across configs A/B in both arms (the mode and
# relaxation keywords change the numerics identically in both arms), so
# the differential is not vacuous. The nondeterministic `Level=5`
# `RealTime` reports are structurally suppressed in both arms via the
# standard `Max Output Level = 4` keyword; `Assembly: (s)` / `Solve:
# (s)` per-iteration timing values (if any) are masked, keeping
# structure and iteration counts compared.
#
# Needs a configured+built serial Elmer tree (ELMER_BUILD — both arms
# share it, so the differential is self-consistent) and the Elmer
# checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:-/home/john/firstmate/data/rs-elmer-t020/scratch/elmer-build}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-filmcompressibility"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T022-A FilmCompressibility proof: work=$WORK build=$ELMER_BUILD"
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
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libfilmcompressibility.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/FilmCompressibility.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the FilmCompressibility module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T filmcompressibility_" "$WORK/port.syms" \
  || { echo "FAIL: filmcompressibility_ not the exported entry"; exit 1; }
grep -qx "T filmcompressibility_init_" "$WORK/port.syms" \
  || { echo "FAIL: filmcompressibility_init_ not exported"; exit 1; }
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
  Simulation Type = Transient
  Timestepping Method = BDF
  BDF Order = 1
  Timestep Sizes(1) = 0.1
  Timestep Intervals(1) = 2
  Output Intervals = 1
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
  Pressure = Real 1.0e6
  Height = Variable Time
    Real MATC "0.5+0.5*tx"
End
Equation 1
  Active Solvers(4) = 1 2 3 4
End
Material 1
  Hydrostatic Pressure = Real 8.0e5
  Thickness = Real 1.0
  Ice Density = Real 900.0
End
Solver 1
  Equation = AllocP
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = Pressure
  Variable DOFs = 1
End
Solver 2
  Equation = AllocH
  Procedure = "AllocateSolver" "AllocateSolver"
  Variable = Height
  Variable DOFs = 1
End
Solver 3
  Equation = FilmEq
  Procedure = "FilmCompressibility" "FilmCompressibility"
  Pressure Variable Name = String "Pressure"
  Height Variable Name = String "Height"
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
  AC Diffusion = Real 1.0e-3
End
Solver 4
  Equation = SaveEq
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = AC
  Operator 1 = max
  Variable 2 = AC
  Operator 2 = min
End
EOF
python3 - "$WORK/caseA.sif" "$WORK/caseB.sif" <<'PYEOF'
import sys
src, dst = sys.argv[1], sys.argv[2]
s = open(src).read()
old = '  AC Diffusion = Real 1.0e-3\n'
new = ('  Differential Mode = Logical True\n'
       '  AC Relaxation Factor = Real 0.5\n')
assert old in s, "caseA diffusion line missing"
s = s.replace(old, new)
open(dst, 'w').write(s)
PYEOF
grep -c "MATC" "$WORK/caseA.sif"
grep -c "Differential Mode" "$WORK/caseB.sif"
grep -c "AC Relaxation Factor" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/FilmCompressibility.so"
cp "$PORT_SO" "$WORK/shadow-port/FilmCompressibility.so"
nm -D --defined-only "$WORK/shadow-port/FilmCompressibility.so" | grep -q " T filmcompressibility_" \
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
norm() { # drop wall-clock lines for the byte-compare; mask (s) timing values
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" | \
    sed -E '/Assembly: \(s\)|Solve: \(s\)/ s/[0-9]+/N/g' || true
}

echo "-- differentials A+B + rerun of A"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-orig/$f" > "$WORK/a0.txt"
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   configs A+B byte-identical (stdout/stderr/scalars) + rerun byte-identical"

echo "-- branch-liveness control (AC Average differs across A/B, both arms)"
for arm in orig port; do
  avga=$(grep "AC Average" "$WORK/runA-$arm/stdout.txt" | head -1)
  avgb=$(grep "AC Average" "$WORK/runB-$arm/stdout.txt" | head -1)
  [ -n "$avga" ] || { echo "FAIL: no AC Average in A($arm)"; exit 1; }
  [ -n "$avgb" ] || { echo "FAIL: no AC Average in B($arm)"; exit 1; }
  [ "$avga" != "$avgb" ] || { echo "FAIL: A($arm) and B($arm) AC Average identical"; exit 1; }
  cmp -s "$WORK/runA-$arm/scalars.dat" "$WORK/runB-$arm/scalars.dat" \
    && { echo "FAIL: A($arm) and B($arm) scalars identical"; exit 1; }
done
echo "   control holds: branches change the numerics in both arms"

echo "FILMCOMPRESSIBILITY PROOF GREEN wall=$((SECONDS - T0))s"
