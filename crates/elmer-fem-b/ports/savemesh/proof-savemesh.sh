#!/bin/bash
# T021-B SaveMesh solver-entry proof: the Rust cdylib replaces SaveMesh.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `savemesh_`);
#  2. two distinct .sif differentials (Save-All-Timesteps config A;
#     single-save + different Mesh Name config B) are byte-identical to
#     the reference except wall-clock lines — stdout, stderr, AND the
#     written mesh trees (`diff -r`);
#  3. a second full process run of config A against the port is
#     byte-identical (idempotent re-entry);
#  4. branch-liveness control: config B takes the other EveryTime branch
#     and a different Mesh Name, so its mesh tree differs from config
#     A's — and both arms agree on it.
#
# Needs a configured+built Elmer tree (ELMER_BUILD, serial config — both
# arms share it, so the differential is self-consistent) and the Elmer
# checkout only for the square.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T021B_LANE:-/home/john/firstmate/data/rs-elmer-t021b}"
WORK="$LANE/scratch/proof-savemesh"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square*
echo "== T021-B SaveMesh proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: gfortran must be 13.3.0"; exit 1;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly (feaadeeac 2026-09-19)"*) ;;
  *) echo "FAIL: rustc +nightly must be 1.100.0-nightly (feaadeeac 2026-09-19)"; exit 1;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="$PORT_DIR/target/release/libsavemesh.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/SaveMesh.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the SaveMesh module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T savemesh_" "$WORK/port.syms" \
  || { echo "FAIL: savemesh_ not the exported entry"; exit 1; }
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
End
Body 1
  Equation = 1
  Material = 1
End
Equation 1
  Active Solvers(1) = 1
End
Material 1
  Density = 1.0
End
Solver 1
  Equation = Save Mesh
  Procedure = "SaveMesh" "SaveMesh"
  Save Mesh = Logical True
  Mesh Name = File "meshA"
  Save Mesh Directory = File "saved/"
  Save All Timesteps = Logical True
End
EOF
sed -e 's/Mesh Name = File "meshA"/Mesh Name = File "meshB"/' \
    -e 's/Save All Timesteps = Logical True/Save All Timesteps = Logical False/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
grep -q 'meshB' "$WORK/caseB.sif"
grep -q 'Logical False' "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/SaveMesh.so"
cp "$PORT_SO" "$WORK/shadow-port/SaveMesh.so"
nm -D --defined-only "$WORK/shadow-port/SaveMesh.so" | grep -q "T savemesh_" \
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
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|ELMER FINISHED AT|STARTED AT" "$1" || true
}
diff_tree() { # <run1> <run2>: byte-compare the written mesh trees
  diff -r "$WORK/$1/saved" "$WORK/$2/saved" > /dev/null 2>&1 || \
    { echo "FAIL: mesh trees $1 vs $2 diverge"; diff -r "$WORK/$1/saved" "$WORK/$2/saved" | head; exit 1; }
}

echo "-- differential A (Save All Timesteps) + rerun (re-entry)"
run_case caseA orig runA-orig
run_case caseA port runA-port
run_case caseA port runA-port2
for f in stdout.txt stderr.txt; do
  norm "$WORK/runA-orig/$f" > "$WORK/a0.txt"
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  diff "$WORK/a0.txt" "$WORK/a1.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config A $f diverges"; diff "$WORK/a0.txt" "$WORK/a1.txt" | head; exit 1; }
done
diff_tree runA-orig runA-port
echo "   config A byte-identical (stdout/stderr/mesh tree)"
for f in stdout.txt stderr.txt; do
  norm "$WORK/runA-port/$f" > "$WORK/a1.txt"
  norm "$WORK/runA-port2/$f" > "$WORK/a2.txt"
  diff "$WORK/a1.txt" "$WORK/a2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: rerun $f diverges"; exit 1; }
done
diff_tree runA-port runA-port2
echo "   rerun byte-identical (re-entry)"

echo "-- differential B (single save, different Mesh Name)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
diff_tree runB-orig runB-port
echo "   config B byte-identical (stdout/stderr/mesh tree)"

echo "-- branch-liveness control (B differs from A, both arms agree)"
diff -r "$WORK/runA-orig/saved" "$WORK/runB-orig/saved" > /dev/null 2>&1 && \
  { echo "FAIL: liveness control dead (A and B mesh trees identical)"; exit 1; }
diff -r "$WORK/runA-port/saved" "$WORK/runB-port/saved" > /dev/null 2>&1 && \
  { echo "FAIL: liveness control dead on port arm"; exit 1; }
echo "   liveness confirmed (mesh trees differ across branches, arms agree)"

echo "SAVEMESH PROOF GREEN ($((SECONDS - T0))s)"
