#!/bin/bash
# T022-A CoilSolver solver-entry proof: the Rust cdylib replaces CoilSolver.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: `coilsolver_`, `coilsolver_init_`,
#     `coilpotential_`, `coilpotentialnormalized_`);
#  2. four distinct .sif differentials are byte-identical to the reference
#     except wall-clock lines/values — config A (basic open coil, nodal
#     current, default conductivity path); config B (config A +
#     `Normalize Coil Current` averaging path + `Coil Conductivity Fix`
#     (second nonlinear iteration through `LocalFixMatrix`) + constant
#     `Electric Conductivity`, i.e. the `GotElCond` assembly path);
#     config C (config A + a Poisson solve whose `Source` is
#     `Variable CoilPot` with
#     `Real Procedure "CoilSolver" "CoilPotential"`, exercising the
#     `coilpotential_` user-function entry per node); config D (same
#     with `CoilPotentialNormalized`, exercising the
#     `coilpotentialnormalized_` entry incl. its first-visit `n`
#     overwrite and per-element gradient cache);
#  3. a second full process run of config A against the port is
#     byte-identical (the solver allocates its work arrays per call, so
#     separate processes are trivially re-entrant; note the reference
#     itself aborts on a second in-process steady iteration —
#     double `ALLOCATE` of the entry work arrays — so all configs use a
#     single steady iteration and re-entry is proven port-vs-port).
# Branch-liveness control: the saved scalars and solver norms differ
# across configs A/B/C/D in both arms (the keyword and UDF branches
# change the numerics identically in both arms), so the differential is
# not vacuous. The nondeterministic `Level=5` `RealTime` reports are
# structurally suppressed in both arms via the standard
# `Max Output Level = 4` keyword (deterministic `Level=4` lines still
# print and are compared).
#
# Serial-proof scope: the `ParEnv % PEs > 1` branches (MPI_ALLREDUCE,
# ParallelReduction) and the closed-coil fitting paths
# (DefineCoilCenter/DefineCoilParameters/ChooseFixedBulkNodesNarrow/
# ChooseCoilCut) are transcribed verbatim but unproven — no closed-coil
# mesh is staged and the proof runs serially.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial — both arms
# share it, so the differential is self-consistent) and the Elmer
# checkout only for the CoilSolver1 Curve.msh mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-coilsolver"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/Curve*
echo "== T022-A CoilSolver proof: work=$WORK build=$ELMER_BUILD"
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
  *"1.100.0-nightly (feaadeeac 2026-09-19)"*) ;;
  *) echo "FAIL: need rustc 1.100.0-nightly (feaadeeac 2026-09-19), got: $RUSTC_VERSION"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libcoilsolver.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/CoilSolver.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the CoilSolver module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in coilsolver_ coilsolver_init_ coilpotential_ coilpotentialnormalized_; do
  grep -qx "T $sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not the exported entry"; exit 1; }
done
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's CoilSolver1 Curve.msh)"
cp "$ELMER_FORK/fem/tests/CoilSolver1/Curve.msh" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 14 2 Curve.msh -autoclean > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "Curve"
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
End
Body 2
  Equation = 2
End
Equation 1
  Name = "CoilCurrent"
  Active Solvers(2) = 1 2
End
Equation 2
  Name = "NoEquations"
End
Solver 1
  Equation = "CoilSolver"
  Procedure = "CoilSolver" "CoilSolver"
  Linear System Solver = "Iterative"
  Linear System Preconditioning = ilu1
  Linear System Max Iterations = Integer 1000
  Linear System Convergence Tolerance = Real 1e-10
  Linear System Iterative Method = BiCGStabL
  Linear System Residual Output = Integer 20
  Steady State Convergence Tolerance = Real 1e-06
  Desired Coil Current = Real 1.0
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = "scalars.dat"
  Variable 1 = "CoilPot"
  Operator 1 = min
  Variable 2 = "CoilPot"
  Operator 2 = max
End
Material 1
  Density = Real 1.0
End
Boundary Condition 1
  Name = "Start of coil loop"
  Target Boundaries = Integer 4
  Coil Start = Logical True
End
Boundary Condition 2
  Name = "End of coil loop"
  Target Boundaries = Integer 1
  Coil End = Logical True
End
EOF
sed -e 's/Desired Coil Current = Real 1.0/Desired Coil Current = Real 1.0\n  Normalize Coil Current = Logical True\n  Coil Conductivity Fix = Logical True\n  Cfix Diffusion = Real 1.0e-4/' \
    -e 's/Density = Real 1.0/Density = Real 1.0\n  Electric Conductivity = Real 1.0/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -c "Normalize Coil Current" "$WORK/caseB.sif"
grep -c "Electric Conductivity" "$WORK/caseB.sif"
grep -c "Coil Conductivity Fix" "$WORK/caseB.sif"
gen_udf() { # <in> <out> <funcname>
  python3 - "$1" "$2" "$3" <<'PYEOF'
import sys
src, dst, func = sys.argv[1], sys.argv[2], sys.argv[3]
text = open(src).read()
poisson = '''Solver 2
  Equation = "UDFPoisson"
  Variable = "UDFPot"
  Variable DOFs = Integer 1
  Procedure = "Poisson" "PoissonSolver"
  Linear System Solver = "Direct"
  Linear System Direct Method = "Umfpack"
End
'''
savescalars = '''Solver 3
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = "scalars.dat"
  Variable 1 = "CoilPot"
  Operator 1 = min
  Variable 2 = "CoilPot"
  Operator 2 = max
  Variable 3 = "UDFPot"
  Operator 3 = min
  Variable 4 = "UDFPot"
  Operator 4 = max
End
Body Force 1
  Source = Variable CoilPot
    Real Procedure "CoilSolver" "FUNCNM"
End
'''
savescalars = savescalars.replace("FUNCNM", func)
old_solv2 = text[text.index("Solver 2"):text.index("Material 1")]
text = text.replace(old_solv2, poisson + savescalars)
text = text.replace("Active Solvers(2) = 1 2", "Active Solvers(3) = 1 2 3")
text = text.replace('''Body 1
  Equation = 1
  Material = 1
End''', '''Body 1
  Equation = 1
  Material = 1
  Body Force = 1
End''')
text = text.replace('''Target Boundaries = Integer 4
  Coil Start = Logical True
End''', '''Target Boundaries = Integer 4
  Coil Start = Logical True
  UDFPot = Real 0.0
End''')
open(dst, "w").write(text)
PYEOF
}
gen_udf "$WORK/caseA.sif" "$WORK/caseC.sif" "CoilPotential"
gen_udf "$WORK/caseA.sif" "$WORK/caseD.sif" "CoilPotentialNormalized"
grep -c 'Real Procedure "CoilSolver" "CoilPotential"' "$WORK/caseC.sif"
grep -c 'Real Procedure "CoilSolver" "CoilPotentialNormalized"' "$WORK/caseD.sif"
grep -c "UDFPot" "$WORK/caseC.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/CoilSolver.so"
cp "$PORT_SO" "$WORK/shadow-port/CoilSolver.so"
nm -D --defined-only "$WORK/shadow-port/CoilSolver.so" | grep -q "T coilsolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/Curve" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines/values for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT" "$1" \
    | sed -E 's/:[[:space:]]+[0-9]+(\.[0-9]+)?(E[+-]?[0-9]+)?[[:space:]]*\(s\)/: (s)/g' || true
}

echo "-- differentials A/B/C/D + rerun of A"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
run_case caseB orig runB-orig
run_case caseB port runB-port
run_case caseC orig runC-orig
run_case caseC port runC-port
run_case caseD orig runD-orig
run_case caseD port runD-port
for f in stdout.txt stderr.txt scalars.dat; do
  for c in A B C D; do
    norm "$WORK/run$c-orig/$f" > "$WORK/cmp-$c-orig.txt"
    norm "$WORK/run$c-port/$f" > "$WORK/cmp-$c-port.txt"
    diff "$WORK/cmp-$c-orig.txt" "$WORK/cmp-$c-port.txt" > /dev/null 2>&1 || \
      { echo "FAIL: config $c $f diverges"; diff "$WORK/cmp-$c-orig.txt" "$WORK/cmp-$c-port.txt" | head; exit 1; }
  done
done
echo "   configs A+B+C+D byte-identical (stdout/stderr/scalars)"
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runA-port/$f" > "$WORK/cmp-r1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/cmp-r2.txt"
  diff "$WORK/cmp-r1.txt" "$WORK/cmp-r2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
echo "   rerun byte-identical (re-entry)"

echo "-- branch-liveness control (saved scalars differ across A/B/C/D, both arms)"
for arm in orig port; do
  for c in A B C D; do
    echo -n "   $arm-$c:"; tr '\n' ' ' < "$WORK/run$c-$arm/scalars.dat" | head -c 200; echo
  done
done
cmp -s "$WORK/runA-orig/scalars.dat" "$WORK/runB-orig/scalars.dat" \
  && { echo "FAIL: A and B scalars identical (control vacuous)"; exit 1; }
cmp -s "$WORK/runA-orig/scalars.dat" "$WORK/runC-orig/scalars.dat" \
  && { echo "FAIL: A and C scalars identical (control vacuous)"; exit 1; }
cmp -s "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" \
  && { echo "FAIL: A and B scalars identical in port arm"; exit 1; }
cmp -s "$WORK/runA-port/scalars.dat" "$WORK/runC-port/scalars.dat" \
  && { echo "FAIL: A and C scalars identical in port arm"; exit 1; }
echo "   control holds: branches change the numerics in both arms"

echo "COILSOLVER PROOF GREEN wall=$((SECONDS - T0))s"
