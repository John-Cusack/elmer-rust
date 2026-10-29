#!/bin/bash
# T022-A GmshOutputReader proof: the Rust cdylib replaces GmshOutputReader.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: only `gmshoutputreader_` + `gmshoutputreader_init_`);
#  2. two distinct .sif differentials that exercise the solver
#     (config A: default search, field inputA.msh; config B: a different
#     ValueList branch `Use Quadrant Tree = Logical True` with 2 steady
#     iterations — i.e. the repeated-invocation run — on different field
#     data inputB.msh) are byte-identical to the reference except
#     wall-clock lines;
#  3. a second full process run of config A is byte-identical
#     (idempotent re-entry through the SAVE state).
#  4. branch-liveness control: config B (different branch + different
#     field data) produces different scalars than config A, and both
#     arms agree on them.
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
WORK="$LANE/scratch/proof-gmshoutputreader"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/square* "$WORK"/input*.msh
echo "== T022-A GmshOutputReader proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: need gfortran 13.3.0, got: $GFORTRAN_VERSION"; exit 1;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly (feaadeeac 2026-09-19)"*) ;;
  *) echo "FAIL: need nightly feaadeeac 2026-09-19, got: $RUSTC_VERSION"; exit 1;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libgmshoutputreader.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- clippy clean"
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly clippy --release -- -D warnings) \
  || { echo "FAIL: clippy warnings"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/GmshOutputReader.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO (build the GmshOutputReader module first)"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
grep -qx "T gmshoutputreader_" "$WORK/port.syms" \
  || { echo "FAIL: gmshoutputreader_ not exported"; exit 1; }
grep -qx "T gmshoutputreader_init_" "$WORK/port.syms" \
  || { echo "FAIL: gmshoutputreader_init_ not exported"; exit 1; }
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's square.grd)"
cp "$ELMER_FORK/fem/tests/PoissonDB/square.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 square > /dev/null)

echo "-- gmsh inputs (clean-room generated from the square mesh nodes/elements)"
python3 - "$WORK" <<'PYEOF'
import sys
W = sys.argv[1]
nodes = []
for line in open(f"{W}/square/mesh.nodes"):
    p = line.split()
    if len(p) < 5:
        continue
    nodes.append((int(p[0]), float(p[2]), float(p[3]), float(p[4])))
elems = []
for line in open(f"{W}/square/mesh.elements"):
    p = line.split()
    if len(p) < 6:
        continue
    assert int(p[2]) == 303, p
    elems.append((int(p[0]), list(map(int, p[3:6]))))
assert len(nodes) > 0 and len(elems) > 0
def write_msh(path, func):
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write("$Nodes\n%d\n" % len(nodes))
        for nid, x, y, z in nodes:
            f.write("%d %.15g %.15g %.15g\n" % (nid, x, y, z))
        f.write("$EndNodes\n")
        f.write("$Elements\n%d\n" % len(elems))
        for eid, ns in elems:
            f.write("%d 2 2 0 1 %d %d %d\n" % (eid, ns[0], ns[1], ns[2]))
        f.write("$EndElements\n")
        f.write("$NodeData\n1\n\"Temperature\"\n1\n0.0\n3\n0\n1\n%d\n" % len(nodes))
        for nid, x, y, z in nodes:
            f.write("%d %.15g\n" % (nid, func(x, y, z)))
        f.write("$EndNodeData\n")
write_msh(f"{W}/inputA.msh", lambda x, y, z: x + 2.0 * y)
write_msh(f"{W}/inputB.msh", lambda x, y, z: 5.0 - 3.0 * x + y)
print("gmsh inputs: %d nodes %d elements" % (len(nodes), len(elems)))
PYEOF

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
  Max Output Level = 4
End
Body 1
  Equation = 1
  Material = 1
End
Body 2
  Equation = 1
  Material = 1
End
Body 3
  Equation = 1
  Material = 1
End
Body 4
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
  Equation = GmshReader
  Procedure = "GmshOutputReader" "GmshOutputReader"
  Variable = Temperature
  Variable DOFs = 1
  Filename = File "inputA.msh"
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = Temperature
  Operator 1 = min
  Variable 2 = Temperature
  Operator 2 = max
End
EOF
sed -e 's/Steady State Max Iterations = 1/Steady State Max Iterations = 2/' \
    -e 's/inputA/inputB/' \
    "$WORK/caseA.sif" > "$WORK/caseB.sif"
python3 - "$WORK/caseB.sif" <<'PYEOF'
import sys
p = sys.argv[1]
s = open(p).read()
s = s.replace('  Filename = File "inputB.msh"\n',
              '  Filename = File "inputB.msh"\n  Use Quadrant Tree = Logical True\n')
open(p, "w").write(s)
PYEOF
grep -q "Use Quadrant Tree = Logical True" "$WORK/caseB.sif"

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/GmshOutputReader.so"
cp "$PORT_SO" "$WORK/shadow-port/GmshOutputReader.so"
nm -D --defined-only "$WORK/shadow-port/GmshOutputReader.so" | grep -q "T gmshoutputreader_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir> <msh>
  local case=$1 arm=$2 run=$3 msh=$4
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/square" "$WORK/$run"/
  cp "$WORK/$msh" "$WORK/$run/square/"
  (cd "$WORK/$run" && \
    ELMER_LIB="$WORK/shadow-$arm" \
    LD_LIBRARY_PATH="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src" \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif > stdout.txt 2> stderr.txt)
  local rc=$?
  echo "$run($arm) exit: $rc"
  [ "$rc" -eq 0 ] || { echo "FAIL: $run($arm) exited $rc"; exit 1; }
}
norm() { # drop wall-clock lines for the byte-compare
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at" "$1" || true
}

echo "-- differential A (default search branch) + rerun (re-entry)"
run_case caseA orig runA-orig inputA.msh
run_case caseA port runA-port inputA.msh
run_case caseA port runA-port2 inputA.msh
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

echo "-- differential B (quadrant-tree branch, 2 iterations, other field data)"
run_case caseB orig runB-orig inputB.msh
run_case caseB port runB-port inputB.msh
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- liveness control (B differs from A on both arms)"
if diff -q "$WORK/runA-orig/scalars.dat" "$WORK/runB-orig/scalars.dat" > /dev/null 2>&1; then
  echo "FAIL: liveness control produced identical scalars"; exit 1
fi
if diff -q "$WORK/runA-port/scalars.dat" "$WORK/runB-port/scalars.dat" > /dev/null 2>&1; then
  echo "FAIL: liveness control produced identical scalars (port)"; exit 1
fi
echo "   liveness control differs on both arms, and both arms agree"

T1=$SECONDS
echo "GMSHOUTPUTREADER PROOF GREEN (wall $((T1-T0))s)"
