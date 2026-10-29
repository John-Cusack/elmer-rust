#!/bin/bash
# T022-A EMWaveSolver proof: the Rust cdylib replaces EMWaveSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: emwavesolver_, emwavesolver_init0_,
#     emwavecalcfields_, emwavecalcfields_init_,
#     emwavecalcfields_init0_);
#  2. three distinct .sif differentials are byte-identical to the
#     reference except wall-clock lines/values:
#       A: 3-D transient shoebox, 2nd-order time + Constant Bulk Matrix,
#          nodal CalcFields (Elfield), driven inlet + damping outlet;
#       B: Time Derivative Order = 1, no constant bulk, 2 nonlinear
#          iterations (re-assembly + DefaultConverged loop), TEM outlet
#          (MATMUL path) — LocalMatrixBC/1st-order-time branch coverage;
#       C: DG/elemental-field path (Init0 solver injection, LocalSol,
#          GlobalSol EL_Var branch) + electric-field time derivatives;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state);
#  4. branch-liveness control: B and C differ from A, and both arms agree
#     on each config.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
# config — both arms share it, so the differential is self-consistent)
# and the Elmer checkout only for the shoebox mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
#
# Divergence note (defined-path rule): the reference fails its own
# iterative solve (IterSolve divergence, exit 1) when the elemental-field
# path (config C shape) is combined with first-order time and no
# constant bulk matrix; the port is not exercised on that
# reference-undefined path. All proof configs below exit 0 on the
# reference.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-emwavesolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/shoebox*
echo "== T022-A EMWaveSolver proof: work=$WORK build=$ELMER_BUILD"
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
  *) echo "FAIL: need rustc nightly 1.100.0, got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libemwavesolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/EMWaveSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the EMWaveSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in emwavesolver_ emwavesolver_init0_ emwavecalcfields_ \
    emwavecalcfields_init_ emwavecalcfields_init0_; do
  grep -qx "T $sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not an exported entry"; exit 1; }
done
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's shoebox_hexas.grd)"
cp "$ELMER_FORK/fem/tests/EMWaveBoxHexas/shoebox_hexas.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 shoebox_hexas > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Check Keywords "Warn"

$ w=2*pi*(2.5e9)
$ eps0=8.854187817e-12
$ mu0=4e-7*pi

Header
  Mesh DB "." "shoebox_hexas"
End

Simulation
  Max Output Level = 4

  Coordinate System = "Cartesian"
  Simulation Type = Transient
  Timestepping Method = Bossak
  Timestep Sizes = $2*pi/(w*50)

  Timestep Intervals = 2
  Steady State Max Iterations = 1

  Output Intervals(1) = 1
End

Constants
  Permittivity of Vacuum = $eps0
  Permeability of Vacuum = $mu0
End

Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End

Material 1
  Relative Permittivity = 1.0
  Relative Permeability = 1.0
End

Equation 1
  Active Solvers(2) = 1 2
End


Solver 1
  Equation = "VectorWave"
  Variable = E

  Optimize Bandwidth = True
  Linear System Symmetric = False
  Procedure = "EMWaveSolver" "EMWaveSolver"

  Quadratic Approximation = Logical False

  Linear System Scaling = True
  Linear System Solver = String "Iterative"
  Linear System Iterative Method = String "bicgstabl"
  BiCGstabl polynomial degree = Integer 4
  Linear System Preconditioning = String "none"
  Linear System ILUT Tolerance = Real 1e-3
  Linear System Max Iterations = 1000
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Abort Not Converged = True
  Linear System Residual Output = 1

  Use Global Mass Matrix = True
  Constant Bulk Matrix = True
End

Solver 2
  Equation = "calcfields"

  Procedure = "EMWaveSolver" "EMWaveCalcFields"
  Linear System Symmetric = False

  Calculate Elemental Fields = Logical False

  Linear System Solver = "Iterative"
  Linear System Preconditioning = None
  Linear System Residual Output = 10
  Linear System Max Iterations = 1000
  Linear System Iterative Method = CG
  Linear System Convergence Tolerance = 1.0e-9

  Constant Bulk Matrix = True

  Exec Solver = before saving
End

Solver 3
  Equation = "save scalars"
  Procedure = "SaveData" "SaveScalars"

  Filename = scalars.dat

  Variable 1 = E
  Operator 1 = min
  Variable 2 = E
  Operator 2 = max
  Variable 3 = Elfield
  Operator 3 = min
  Variable 4 = Elfield
  Operator 4 = max

  Exec Solver = after saving
End

Body Force 1
End


Boundary Condition 1
  Name = "Inlet"
  Target Boundaries(1) = 5

  E {e} 1 = Real 0
  E {e} 2 = Variable "time"
    Real MATC "sin(w*tx)"
  E {e} 3 = Real 0
End

Boundary Condition 2
  Name = "Outlet"
  Target Boundaries(1) = 6

  Electric Damping Coefficient = Real $ sqrt(eps0 * mu0)
End

Boundary Condition 3
  Name = "Sides"
  Target Boundaries(4) = 1 2 3 4

  E {e} 1 = Real 0
  E {e} 2 = Real 0
  E {e} 3 = Real 0
End
EOF
# Config B: first-order time, no constant bulk, 2 nonlinear iterations,
# TEM outlet (GotTem MATMUL path).
sed -e 's/  Quadratic Approximation = Logical False/  Quadratic Approximation = Logical False\n  Time Derivative Order = 1\n  Nonlinear System Max Iterations = Integer 2/' \
    -e 's/  Constant Bulk Matrix = True/  Constant Bulk Matrix = Logical False/' \
    -e 's|  Electric Damping Coefficient = Real $ sqrt(eps0 \* mu0)|  TEM Potential = Real 1.0\n  Electric Damping Coefficient = Real $ sqrt(eps0 * mu0)|' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
# Config C: DG/elemental-field path (Init0 solver injection) + electric
# field time derivatives; Elfield E and dEdt saved for observability.
python3 - "$WORK/caseA.sif" "$WORK/caseC.sif" <<'PYEOF'
import sys
src, dst = sys.argv[1], sys.argv[2]
t = open(src).read()
t = t.replace("  Calculate Elemental Fields = Logical False\n",
              "  Calculate Electric field derivatives = Logical True\n")
t = t.replace("""  Variable 4 = Elfield
  Operator 4 = max""",
              """  Variable 4 = Elfield
  Operator 4 = max
  Variable 5 = Elfield E
  Operator 5 = min
  Variable 6 = dEdt
  Operator 6 = max""")
open(dst, "w").write(t)
PYEOF
grep -c "Time Derivative Order = 1" "$WORK/caseB.sif"
grep -c "TEM Potential" "$WORK/caseB.sif"
grep -c "Calculate Electric field derivatives" "$WORK/caseC.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/EMWaveSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/EMWaveSolver.so"
nm -D --defined-only "$WORK/shadow-port/EMWaveSolver.so" | grep -q "T emwavesolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/shoebox_hexas" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines; mask per-iteration timing values only
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT" "$1" \
    | sed -E 's/[0-9]+\.[0-9]+(E[+-]?[0-9]+)? \(s\)/(s)/g' || true
}

echo "-- differential A (2nd-order time + constant bulk + nodal fields) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stderr.txt scalars.dat; do
  diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/runA-orig/$f" "$WORK/runA-port/$f" | head; exit 1; }
done
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a1.txt"
diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
  { echo "FAIL: config A stdout diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
echo "   config A byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (SAVE re-entry)"

echo "-- differential B (1st-order time, no bulk, TEM outlet)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- differential C (DG/elemental fields + time derivatives)"
run_case caseC orig runC-orig
run_case caseC port runC-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runC-orig/$f" > "$WORK/c1.txt"
  norm "$WORK/runC-port/$f" > "$WORK/c2.txt"
  diff "$WORK/c1.txt" "$WORK/c2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config C $f diverges"; diff "$WORK/c1.txt" "$WORK/c2.txt" | head; exit 1; }
done
echo "   config C byte-identical (stdout/stderr/scalars)"

echo "-- branch-liveness control (variants differ, arms agree)"
norm "$WORK/runA-port/stdout.txt" > "$WORK/a-live.txt"
norm "$WORK/runB-port/stdout.txt" > "$WORK/b2.txt"
norm "$WORK/runC-port/stdout.txt" > "$WORK/c2.txt"
diff "$WORK/a-live.txt" "$WORK/b2.txt" > /dev/null 2>&1 && \
  { echo "FAIL: config B output identical to A (branch not live)"; exit 1; }
diff "$WORK/a-live.txt" "$WORK/c2.txt" > /dev/null 2>&1 && \
  { echo "FAIL: config C output identical to A (branch not live)"; exit 1; }
diff "$WORK/runA-orig/scalars.dat" "$WORK/runB-orig/scalars.dat" > /dev/null 2>&1 && \
  { echo "FAIL: config B scalars identical to A (branch not live)"; exit 1; }
echo "   variants differ on both arms"

echo "EMWAVESOLVER PROOF GREEN (wall $((SECONDS - T0))s)"
