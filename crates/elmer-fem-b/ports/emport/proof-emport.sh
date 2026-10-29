#!/bin/bash
# T022-A EMPort proof: the Rust cdylib replaces EMPort.so.
#
# Proves, on pinned gfortran 13.3.0 + campaign nightly:
#  1. the cdylib exports exactly the reference .so's dynamic symbols
#     (`nm -D` equality: emportsolver_, emportsolver_init0_,
#     emportpotential_, emportpotential_init_);
  #  2. two distinct .sif differentials are byte-identical to the reference
  #     except wall-clock lines/values (config A: base eigen formulation;
  #     config B: potentials formulation + nodal post + Eigenfunction Index 2);
  #  3. a second full process run of config A is byte-identical
  #     (SAVE re-entry);
  #  4. branch-liveness control: configs A and B produce different output,
  #     and both arms agree on each.
  #
  # Branch notes: config B exercises LocalMatrix_With_Potentials, EMPortPost,
  # the Eigenfunction Index 2 path, and the Init0 Port Ground presets. The
  # MaxPort >= 1 per-element port filter and the MaxPort > 1 multi-matrix
  # surgery need port-tagged BCs / a multi-port 3-D mesh (upstream
  # EM_port_eigen_3D); they are transcribed statement-for-statement but not
  # exercised here. The ParEnv % PEs > 1 paths are guarded out exactly as in
  # the original (serial proof). The reference run is deterministic serial
  # ARPACK/UMFPACK output; no defined-path divergence was observed.
  #
  # Needs a configured+built Elmer tree (ELMER_BUILD, serial or campaign
  # config — both arms share it, so the differential is self-consistent)
  # and the Elmer checkout only for the port.grd mesh input.
# Honest halts: any divergence is a FAIL.
# NEVER edit this script to make it pass — fix the port or the wrappers.
set -euo pipefail
PORT_DIR="$(cd "$(dirname "$0")" && pwd)"
ELMER_BUILD="${ELMER_BUILD:?set ELMER_BUILD to a built Elmer tree}"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
LANE="${RS_ELMER_T022A_LANE:-/home/john/firstmate/data/rs-elmer-t022a}"
WORK="$LANE/scratch/proof-emport"
mkdir -p "$WORK"
rm -rf "$WORK"/run* "$WORK"/shadow-* "$WORK"/case*.sif "$WORK"/port*
echo "== T022-A EMPort proof: work=$WORK build=$ELMER_BUILD"
T0=$SECONDS

echo "-- pin toolchain"
GFORTRAN_VERSION="$(gfortran --version | head -1)"
echo "   $GFORTRAN_VERSION"
case "$GFORTRAN_VERSION" in
  *"13.3.0"*) ;;
  *) echo "FAIL: gfortran pin mismatch"; exit 1 ;;
esac
RUSTC_VERSION="$(rustc +nightly --version)"
echo "   $RUSTC_VERSION"
case "$RUSTC_VERSION" in
  *"1.100.0-nightly"*) ;;
  *) echo "FAIL: rustc pin mismatch"; exit 1 ;;
esac

echo "-- build port cdylib"
test -f "$ELMER_BUILD/fmodules/types.mod" || { echo "FAIL: no .mod files under ELMER_BUILD/fmodules"; exit 1; }
(cd "$PORT_DIR" && ELMER_MODDIR="$ELMER_BUILD/fmodules" cargo +nightly build --release -q)
PORT_SO="${CARGO_TARGET_DIR:-$PORT_DIR/target}/release/libemport.so"
test -f "$PORT_SO" || { echo "FAIL: no $PORT_SO"; exit 1; }

echo "-- symbol surface equality (nm -D)"
REF_SO="$ELMER_BUILD/fem/src/modules/EMPort.so"
test -f "$REF_SO" || { echo "FAIL: no $REF_SO"; exit 1; }
nm -D --defined-only "$REF_SO" | awk '{print $2, $3}' | sort > "$WORK/ref.syms"
nm -D --defined-only "$PORT_SO" | awk '{print $2, $3}' | sort > "$WORK/port.syms"
cat "$WORK/ref.syms"
diff "$WORK/ref.syms" "$WORK/port.syms" \
  || { echo "FAIL: dynamic symbol surface differs"; exit 1; }
for sym in emportsolver_ emportsolver_init0_ emportpotential_ emportpotential_init_; do
  grep -qx "T $sym" "$WORK/port.syms" \
    || { echo "FAIL: $sym not an exported entry"; exit 1; }
done
echo "   symbol surfaces match exactly"

echo "-- mesh (ElmerGrid from the Elmer checkout's EM_port_eigen/port.grd)"
cp "$ELMER_FORK/fem/tests/EM_port_eigen/port.grd" "$WORK/"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 port > /dev/null)

echo "-- configs (clean-room .sif files below)"
cat > "$WORK/caseA.sif" <<'EOF'
Header
  Mesh DB "." "port"
End
Simulation
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Max Output Level = 4
End
Constants
  Permeability of Vacuum = 1.0
  Permittivity of Vacuum = 1.0
End
Body 1
  Equation = 1
  Material = 1
End
Body 2
  Equation = 1
  Material = 2
End
Material 1
  Name = "higher_eps"
  Relative Reluctivity = 1.0
  Relative Permittivity = 2.25
End
Material 2
  Name = "Vacuum"
  Relative Reluctivity = 1.0
  Relative Permittivity = 1.0
End
Equation 1
  Active Solvers(1) = 1
End
Solver 1
  Equation = "Port mode"
  Procedure = "EMPort" "EMPortSolver"
  Variable Output = False
  Angular Frequency = Real 3.0
  Use Piola Transform = True
  Linear System Solver = "Direct"
  Linear System Direct Method = Umfpack
  Eigen System Shift Automatic = Logical True
  Eigen System Sorting = String "smallest real part"
  Eigen Analysis = True
  Eigen System Values = 30
  Eigen System Convergence Tolerance = Real 1.0e-10
  Eigen System Compute Residuals = True
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Save EigenValues = True
  Show Norm Index = 1
End
Boundary Condition 1
  Target Boundaries(6) = 1 2 3 4 5 6
  Eport re = Real 0
  Eport im = Real 0
  Eport re {e} = Real 0.0
  Eport im {e} = Real 0.0
End
EOF
cat > "$WORK/caseB.sif" <<'EOF'
Header
  Mesh DB "." "port"
End
Simulation
  Coordinate System = Cartesian
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Max Output Level = 4
End
Constants
  Permeability of Vacuum = 1.0
  Permittivity of Vacuum = 1.0
End
Body 1
  Equation = 1
  Material = 1
End
Body 2
  Equation = 1
  Material = 2
End
Material 1
  Name = "higher_eps"
  Relative Reluctivity = 1.0
  Relative Permittivity = 2.25
End
Material 2
  Name = "Vacuum"
  Relative Reluctivity = 1.0
  Relative Permittivity = 1.0
End
Equation 1
  Active Solvers(1) = 1
End
Solver 1
  Equation = "Port mode"
  Procedure = "EMPort" "EMPortSolver"
  Variable Output = False
  Angular Frequency = Real 3.0
  Use Piola Transform = True
  Use Potential = Logical True
  Linear System Solver = "Direct"
  Linear System Direct Method = Umfpack
  Eigen System Shift Automatic = Logical True
  Eigen System Sorting = String "smallest real part"
  Calculate Nodal Field = Logical True
  Eigenfunction Index = Integer 2
  Eigen Analysis = True
  Eigen System Values = 30
  Eigen System Convergence Tolerance = Real 1.0e-10
  Eigen System Compute Residuals = True
  post: Linear System Solver = Iterative
  post: Linear System Convergence Tolerance = real 1e-8
  post: Linear System Iterative Method = CG
  post: Linear System Residual Output = integer 10
  post: Linear System Max Iterations = integer 2000
  post: Linear System Preconditioning = ILU0
End
Solver 2
  Equation = "SaveScalars"
  Procedure = "SaveData" "SaveScalars"
  Filename = scalars.dat
  Save EigenValues = True
  Show Norm Index = 1
End
Boundary Condition 1
  Target Boundaries(6) = 1 2 3 4 5 6
  Port Ground = True
End
EOF

echo "-- shadow module dirs (port .so under its own name, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"
  mkdir -p "$WORK/shadow-$arm"
  ln -s "$ELMER_BUILD"/fem/src/modules/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-port/EMPort.so"
cp "$PORT_SO" "$WORK/shadow-port/EMPort.so"
nm -D --defined-only "$WORK/shadow-port/EMPort.so" | grep -q "T emportsolver_" \
  || { echo "FAIL: staged port .so lost its entry"; exit 1; }

run_case() { # <case> <arm> <rundir>
  local case=$1 arm=$2 run=$3
  rm -rf "$WORK/$run"
  mkdir -p "$WORK/$run"
  cp "$WORK/$case.sif" "$WORK/$run/case.sif"
  cp -r "$WORK/port" "$WORK/$run"/
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
    | sed -E 's/: +[0-9]+\.[0-9]+ *\(s\)/: (s)/g' || true
}

echo "-- differential A (base eigen formulation) + rerun (re-entry)"
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

echo "-- differential B (potentials + nodal post + port filter)"
run_case caseB orig runB-orig
run_case caseB port runB-port
for f in stdout.txt stderr.txt scalars.dat; do
  norm "$WORK/runB-orig/$f" > "$WORK/b1.txt"
  norm "$WORK/runB-port/$f" > "$WORK/b2.txt"
  diff "$WORK/b1.txt" "$WORK/b2.txt" > /dev/null 2>&1 || \
    { echo "FAIL: config B $f diverges"; diff "$WORK/b1.txt" "$WORK/b2.txt" | head; exit 1; }
done
echo "   config B byte-identical (stdout/stderr/scalars)"

echo "-- liveness control (A and B must differ; both arms must agree)"
norm "$WORK/runA-orig/stdout.txt" > "$WORK/a0.txt"
norm "$WORK/runB-orig/stdout.txt" > "$WORK/b0.txt"
diff -q "$WORK/a0.txt" "$WORK/b0.txt" > /dev/null 2>&1 \
  && { echo "FAIL: liveness control did not diverge"; exit 1; }
echo "   A/B outputs differ on both arms as expected"

echo "EMPORT PROOF GREEN"
echo "elapsed: $((SECONDS - T0))s"
