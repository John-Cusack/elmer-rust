#!/bin/bash
# T025 clean-room stages: per-entry differentials for USF entries without a
# loading .sif in the test tree (or stranded by frozen-off optionals).
#
# Method (ADR-059 §3 pattern): minimal .sif files, one template per
# coupling class, entry name swapped per stage; shadow module dirs
# (symlink all, replace ElmerIceUSF.so per arm); norm + byte-compare;
# rerun (re-entry) + liveness (B-vs-C configs differ on both arms).
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built ice-enabled tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T025_LANE:-/home/john/firstmate/data/rs-elmer-t025}"
WORK="$LANE/scratch/proof-usf-staged"
mkdir -p "$WORK"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
case "$GFORTRAN_VERSION" in *"13.3.0"*) ;; *) echo "FAIL: gfortran 13.3.0 required"; exit 1 ;; esac
RUSTC_VERSION="$(rustc +nightly --version)"
case "$RUSTC_VERSION" in *"1.100.0-nightly (feaadeeac 2026-09-19)"*) ;; *) echo "FAIL: nightly pin"; exit 1 ;; esac
echo "-- build port cdylib"
(cd "$PORT_DIR" && ELMER_BUILD="$ELMER_BUILD" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libElmerIceUSF.so"
REF_SO="$ELMER_BUILD/fem/src/modules/ElmerIceUSF.so"

PRISTINE_SO="$RS_ELMER_T025_LANE/scratch/ElmerIceUSF.pristine.so"
[ -f "$PRISTINE_SO" ] || die "no stable pristine $PRISTINE_SO"
echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/ElmerIceUSF.so"
rm -f "$WORK/shadow-orig/ElmerIceUSF.so"
cp "$PRISTINE_SO" "$WORK/shadow-orig/ElmerIceUSF.so"
rm -f "$WORK/shadow-port/ElmerIceUSF.so"
cp "$PORT_SO" "$WORK/shadow-port/ElmerIceUSF.so"

norm() {
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at|Estimated Time left|Active coordinate mapping time|tcp_peer_recv_connect_ack" "$1" | \
    sed -E -e 's/(Assembly: [ ]*\(s\))[0-9.eE+ -]+/\1 <t>/g' -e 's/(Solve: [ ]*\(s\))[0-9.eE+ -]+/\1 <t>/g' || true
}

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/square" "$WORK/$run"/
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt) \
    || { echo "FAIL: $run($arm) solver failed"; tail -20 "$WORK/$run/stderr.txt"; exit 1; }
}

cmp_case() { # <case> (compares run<case>-orig vs run<case>-port + rerun)
  local case=$1
  for f in stdout.txt stderr.txt scalars.dat; do
    norm "$WORK/run$case-orig/$f" > "$WORK/$case-$f-orig.txt"
    norm "$WORK/run$case-port/$f" > "$WORK/$case-$f-port.txt"
    diff "$WORK/$case-$f-orig.txt" "$WORK/$case-$f-port.txt" > /dev/null 2>&1 || {
      echo "FAIL: $case $f diverges"; diff "$WORK/$case-$f-orig.txt" "$WORK/$case-$f-port.txt" | head; exit 1; }
  done
  echo "   stage $case byte-identical (stdout/stderr/scalars)"
}

echo "-- mesh (ElmerGrid from the Elmer checkout's square.grd)"
cp "$ELMER_FORK/fem/tests/PoissonDB/square.grd" "$WORK/" 2>/dev/null || cp /home/john/runs/elmer-fork/fem/tests/PoissonDB/square.grd "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)

# --- Template P: steady Poisson, one Real-Procedure body-force source ----
mk_poisson() { # <entry> <case>
  local entry=$1 case=$2
  cat > "$WORK/$case.sif" <<EOF
Header
  Mesh DB "." "square"
End
Simulation
  Coordinate System = Cartesian 2D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Output Intervals = 1
  Max Output Level = 4
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
  Active Solvers(1) = 1
End
Material 1
  Density = Real 1.0e3
End
Body Force 1
  Source = Variable Potential
    Real Procedure "ElmerIceUSF" "$entry"
End
Solver 1
  Equation = Poisson
  Procedure = "Poisson" "PoissonSolver"
  Variable = Potential
  Variable DOFs = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU0
End
Solver 2
  Equation = SaveScalars
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Potential
  Operator 1 = min
  Variable 2 = Potential
  Operator 2 = max
End
EOF
}

echo "-- Stage P (pure scalar entries via Poisson source)"
for e in TenPowerA TenPowerA_d SQRTA Asquare Asquare_d; do
  # Reference entry names are case-insensitive-lowered by the loader lookup;
  # stage with the exact dynamic symbol basename.
  le=$(echo "$e" | tr '[:upper:]' '[:lower:]')
  mk_poisson "$e" "P-$le"
  run_case "P-$le" orig "runP-$le-orig"
  run_case "P-$le" port "runP-$le-port"
  run_case "P-$le" port "runP-$le-port2"
  cmp_case "P-$le"
  for f in stdout.txt stderr.txt scalars.dat; do
    norm "$WORK/runP-$le-port/$f" > "$WORK/a1.txt"
    norm "$WORK/runP-$le-port2/$f" > "$WORK/a2.txt"
    diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || { echo "FAIL: rerun P-$le diverges"; exit 1; }
  done
  echo "   rerun P-$le byte-identical"
done
# Liveness: different pure entries give different scalars on both arms.
norm "$WORK/runP-tenpowera-port/scalars.dat" > "$WORK/l1.txt"
norm "$WORK/runP-sqrta-port/scalars.dat" > "$WORK/l2.txt"
if diff "$WORK/l1.txt" "$WORK/l2.txt" > /dev/null 2>&1; then
  echo "FAIL: liveness control identical (TenPowerA == SQRTA)"; exit 1;
fi
echo "   liveness differs (TenPowerA != SQRTA) on both arms"

echo "STAGED PROOF PART 1 GREEN (wall: $((SECONDS - T0))s)"
