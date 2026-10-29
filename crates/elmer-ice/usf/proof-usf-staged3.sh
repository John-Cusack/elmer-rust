#!/bin/bash
# proof-usf-staged3.sh -- Stages S2 (BC/bulk-keyed entries) and S3 (stress-needing
# entries via a ComputeDevStress clone).
#
# Method (ADR-059 S3 pattern): clean-room template per entry class, NO loading
# .sif in the test tree. Lane work only; never /tmp. Requires ELMER_BUILD
# (ice-enabled build) and RS_ELMER_T025_LANE (lane dir).
set -euo pipefail
: "${ELMER_BUILD:?set ELMER_BUILD to the ice-enabled build dir}"
: "${RS_ELMER_T025_LANE:?set RS_ELMER_T025_LANE to the lane dir}"
export HWLOC_COMPONENTS=-gl
die(){ echo "FAIL: $*" >&2; exit 1; }
[ -x "$ELMER_BUILD/fem/src/ElmerSolver" ] || die "no ElmerSolver under $ELMER_BUILD"
ELMER_FORK="${ELMER_FORK:-/home/john/runs/elmer-fork}"
REF_MOD="$ELMER_BUILD/fem/src/modules"
REF_SO="$REF_MOD/ElmerIceUSF.so"
[ -f "$REF_SO" ] || die "no reference $REF_SO"
PORT_SO="$RS_ELMER_T025_LANE/elmer-rust/crates/elmer-ice/usf/target/release/libElmerIceUSF.so"
[ -f "$PORT_SO" ] || die "no port $PORT_SO"
SOLVER_LIBS="$ELMER_BUILD/fem/src:$ELMER_BUILD/matc/src:$ELMER_BUILD/fhutiter/src:$ELMER_BUILD/umfpack/src"

WORK="$RS_ELMER_T025_LANE/scratch/proof-usf-staged3"
rm -rf "$WORK"; mkdir -p "$WORK"

PRISTINE_SO="$RS_ELMER_T025_LANE/scratch/ElmerIceUSF.pristine.so"
[ -f "$PRISTINE_SO" ] || die "no stable pristine $PRISTINE_SO"
echo "-- shadow module dirs (port replaces, rest symlinked)"
for arm in orig port; do
  rm -rf "$WORK/shadow-$arm"; mkdir -p "$WORK/shadow-$arm"
  ln -s "$REF_MOD"/*.so "$WORK/shadow-$arm"/
  ln -s "$ELMER_BUILD"/fem/src/elements.def "$WORK/shadow-$arm"/
  ln -s "$ELMER_FORK"/fem/src/SOLVER.KEYWORDS "$WORK/shadow-$arm"/
done
rm -f "$WORK/shadow-orig/ElmerIceUSF.so"
cp "$PRISTINE_SO" "$WORK/shadow-orig/ElmerIceUSF.so"
rm -f "$WORK/shadow-port/ElmerIceUSF.so"
cp "$PORT_SO" "$WORK/shadow-port/ElmerIceUSF.so"

norm(){
  grep -v -E "Elapsed REAL time|SOLVER TOTAL TIME|FINISHED AT|STARTED AT|File started at|Estimated [Tt]ime left|Active coordinate mapping time|tcp_peer_recv_connect_ack" "$1" \
    | sed -E -e 's/(Assembly: [ ]*\(s\))[0-9.eE+ -]+/\1 <t>/g' -e 's/(Solve: [ ]*\(s\))[0-9.eE+ -]+/\1 <t>/g' || true
}

TISMIP=/home/john/runs/elmer-work/orig/elmerice/Tests/Friction_Weertman
TDEV=/home/john/runs/elmer-work/orig/elmerice/Tests/ComputeDevStress
[ -f "$TISMIP/rectangle.grd" ] || die "no ISMIP rectangle.grd"
cp "$TISMIP/rectangle.grd" "$WORK/rectangle.grd"
(cd "$WORK" && "$ELMER_BUILD/elmergrid/src/ElmerGrid" 1 2 rectangle.grd >/dev/null 2>&1) \
  || die "ElmerGrid failed"

run_pair(){ # <case>: orig+port+rereun, diff stdout/stderr/scalars.dat
  local case="$1"
  for arm in orig port; do
    local r="$WORK/run$case-$arm"
    rm -rf "$r"; mkdir -p "$r"
    cp "$WORK/$case.sif" "$r/case.sif"
    cp -r "$WORK/rectangle" "$r/rectangle"
    ( cd "$r" && ELMER_LIB="$WORK/shadow-$arm" LD_LIBRARY_PATH="$SOLVER_LIBS" timeout 570 \
      "$ELMER_BUILD/fem/src/ElmerSolver" case.sif >stdout.txt 2>stderr.txt ) \
      || { echo "FAIL: $case/$arm solver failed"; tail -n 5 "$r/stderr.txt"; exit 1; }
  done
  norm "$WORK/run$case-orig/stdout.txt" > "$WORK/$case-stdout-orig.txt"
  norm "$WORK/run$case-port/stdout.txt" > "$WORK/$case-stdout-port.txt"
  diff "$WORK/$case-stdout-orig.txt" "$WORK/$case-stdout-port.txt" > /dev/null 2>&1 \
    || { echo "FAIL: $case stdout diverges"; diff "$WORK/$case-stdout-orig.txt" "$WORK/$case-stdout-port.txt" | head -n 10; exit 1; }
  norm "$WORK/run$case-orig/stderr.txt" > "$WORK/$case-stderr-orig.txt"
  norm "$WORK/run$case-port/stderr.txt" > "$WORK/$case-stderr-port.txt"
  diff "$WORK/$case-stderr-orig.txt" "$WORK/$case-stderr-port.txt" > /dev/null 2>&1 \
    || { echo "FAIL: $case stderr diverges"; diff "$WORK/$case-stderr-orig.txt" "$WORK/$case-stderr-port.txt" | head -n 10; exit 1; }
  diff "$WORK/run$case-orig/scalars.dat" "$WORK/run$case-port/scalars.dat" > /dev/null 2>&1 \
    || { echo "FAIL: $case scalars.dat diverges"; diff "$WORK/run$case-orig/scalars.dat" "$WORK/run$case-port/scalars.dat" | head -n 10; exit 1; }
  ( cd "$WORK/run$case-port" && ELMER_LIB="$WORK/shadow-port" LD_LIBRARY_PATH="$SOLVER_LIBS" timeout 570 \
    "$ELMER_BUILD/fem/src/ElmerSolver" case.sif >stdout2.txt 2>&1 ) \
    || { echo "FAIL: $case port rerun failed"; exit 1; }
  norm "$WORK/run$case-port/stdout2.txt" > "$WORK/$case-stdout-port2.txt"
  diff "$WORK/$case-stdout-orig.txt" "$WORK/$case-stdout-port2.txt" > /dev/null 2>&1 \
    || { echo "FAIL: $case port rerun stdout diverges"; exit 1; }
  echo "   stage $case byte-identical (stdout/stderr/scalars.dat)"
  echo "   rerun $case byte-identical"
}

mk_base(){ # <case> <entry>: ISMIP clone with Slip line swapped + SaveScalars
  local case="$1" entry="$2"
  sed -e "s/Sliding_Weertman/$entry/" "$TISMIP/ismip_weertman.sif" > "$WORK/$case.sif"
  sed -i 's/Active Solvers(3) = 1 2 3/Active Solvers(4) = 1 2 3 4/' "$WORK/$case.sif"
  cat >> "$WORK/$case.sif" <<'EOF'

Solver 4
  Equation = SaveScalars
  Procedure = File "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = dSdt
  Operator 1 = min
  Variable 2 = dSdt
  Operator 2 = max
End
EOF
}

echo "-- Stage S2 (BC/bulk-keyed entries via ISMIP clone)"
# -- FreeSlipShelves: mask + formulation keys (dSdt slot keeps BetaIn small)
mk_base S-freeslipshelves FreeSlipShelves
sed -i 's|Slip Coefficient 2 =  Variable Coordinate 1|Slip Coefficient 2 =  Variable dSdt|' "$WORK/S-freeslipshelves.sif"
python3 - "$WORK/S-freeslipshelves.sif" <<'PYEOF'
import sys
p=sys.argv[1]; s=open(p).read()
s=s.replace("Boundary Condition 1\n  Target Boundaries = 1","Boundary Condition 1\n  Target Boundaries = 1\n  FreeSlipShelves mask name = String \"dSdt\"\n  FreeSlipShelves beta formulation = String \"Power\"",1)
open(p,'w').write(s)
PYEOF
run_pair S-freeslipshelves
echo "STAGE-S2A-GREEN $(date -u +%FT%TZ)"
# -- Calculate_HAF: thickness var name + rho ratio
mk_base S-calculate_haf Calculate_HAF
python3 - "$WORK/S-calculate_haf.sif" <<'PYEOF'
import sys
p=sys.argv[1]; s=open(p).read()
s=s.replace("Boundary Condition 1\n  Target Boundaries = 1","Boundary Condition 1\n  Target Boundaries = 1\n  Haf Thickness = String \"dSdt\"\n  Haf rhoRatio = Real 0.913\n  Haf GroundedMask = String \"dSdt\"",1)
open(p,'w').write(s)
PYEOF
run_pair S-calculate_haf
# -- caffeGetViscosity via its natural Viscosity slot (bulk context)
mk_base S-caffegetviscosity caffeGetViscosity
python3 - "$WORK/S-caffegetviscosity.sif" <<'PYEOF'
import sys
p=sys.argv[1]; s=open(p).read()
s=s.replace('Viscosity Model = String "power law"\n  Viscosity = Real $eta','Viscosity = Variable Coordinate 1\n    Real Procedure "ElmerIceUSF" "caffeGetViscosity"',1)
s=s.replace('Real Procedure "ElmerIceUSF" "caffeGetViscosity"\n  Weertman','Real Procedure "ElmerIceUSF" "Sliding_Weertman"\n  Weertman',1)
s=s.replace('Critical Shear Rate = Real 1.0e-10','Critical Shear Rate = Real 1.0e-10\n  Activation Energies(2) = Real 6.0e4 6.0e4\n  Arrhenius Factors(2) = Real 1.0 1.0\n  Limit Temperature = Real -10.0\n  Anisotropic Enhancement Factor = Real 1.0\n  Critical Enhancement factor = Real 1.0e-4\n  Critical Strain Rate = Real 1.0e-15\n  Temperature Name = String "CaffeTemp"',1)
open(p,'w').write(s)
PYEOF
run_pair S-caffegetviscosity
echo "STAGE-S2B-GREEN $(date -u +%FT%TZ)"
# -- RelativeTemperature / ArrheniusFactor via bulk Flow BodyForce slot
for e in RelativeTemperature ArrheniusFactor; do
  le=$(echo "$e" | tr '[:upper:]' '[:lower:]')
  sed -e "s/Sliding_Weertman/$e/" "$TISMIP/ismip_weertman.sif" > "$WORK/S-$le.sif"
  sed -i 's/Active Solvers(3) = 1 2 3/Active Solvers(4) = 1 2 3 4/' "$WORK/S-$le.sif"
  cat >> "$WORK/S-$le.sif" <<'EOF'

Solver 4
  Equation = SaveScalars
  Procedure = File "SaveData" "SaveScalars"
  Filename = scalars.dat
  Variable 1 = dSdt
  Operator 1 = min
  Variable 2 = dSdt
  Operator 2 = max
End
EOF
  python3 - "$WORK/S-$le.sif" "$e" <<'PYEOF'
import sys
p,n=sys.argv[1],sys.argv[2]; s=open(p).read()
s=s.replace('Flow BodyForce 1 = Real 0.0','Flow BodyForce 1 = Variable "Normal Vector"\n    Real Procedure "ElmerIceUSF" "%s"'%n,1)
s=s.replace('! No constant Needed','Gas Constant = Real 8.314\n  Clausius Clapeyron Constant = Real 9.8e-08',1)
s=s.replace('Critical Shear Rate = Real 1.0e-10','Critical Shear Rate = Real 1.0e-10\n  Glen Exponent = Real 3.0\n  Limit Temperature = Real -10.0\n  Rate Factor 1 = Real 1.0\n  Rate Factor 2 = Real 1.0\n  Scale Rate Factors = Logical False\n  Activation Energy 1 = Real 6.0e4\n  Activation Energy 2 = Real 6.0e4\n  Glen Enhancement Factor = Real 1.0',1)
open(p,'w').write(s)
PYEOF
  run_pair "S-$le"
done
echo "STAGE-S2C-GREEN $(date -u +%FT%TZ)"
# -- SourceCalc: melt flags off + dummy hydraulic-potential solver
mk_base S-sourcecalc SourceCalc
python3 - "$WORK/S-sourcecalc.sif" <<'PYEOF'
import sys
p=sys.argv[1]; s=open(p).read()
s=s.replace('Flow BodyForce 1 = Real 0.0','Flow BodyForce 1 = Variable Coordinate 1\n    Real Procedure "ElmerIceUSF" "SourceCalc"\n  Internal Melt = Logical False\n  Surface Melt = Logical False',1)
s=s.replace('Real Procedure "ElmerIceUSF" "SourceCalc"\n  Weertman','Real Procedure "ElmerIceUSF" "Sliding_Weertman"\n  Weertman',1)
s=s.replace('Active Solvers(4) = 1 2 3 4','Active Solvers(5) = 1 2 3 4 5',1)
s=s.replace('Boundary Condition 1\n  Target Boundaries = 1','Boundary Condition 1\n  Target Boundaries = 1\n  hydraulic potential = Real 1.0',1)
open(p,'w').write(s)
PYEOF
cat >> "$WORK/S-sourcecalc.sif" <<'EOF'

Solver 5
  Equation = "Hydraulic"
  Procedure = "Poisson" "PoissonSolver"
  Variable = "hydraulic potential"
  Variable DOFs = 1
  Linear System Solver = Direct
  Linear System Direct Method = umfpack
  Steady State Convergence Tolerance = Real 1.0e-5
End
EOF
run_pair S-sourcecalc
echo "STAGE-S2D-GREEN $(date -u +%FT%TZ)"
echo "-- Stage S3 (stress-needing entries via ComputeDevStress clone)"
mk_stress(){ # <case> <entry>: DevStress clone + Normal solver + ISMIP BC1
  local case="$1" entry="$2"
  cp "$TDEV/computeDevStress.sif" "$WORK/$case.sif"
  python3 - "$WORK/$case.sif" "$entry" <<'PYEOF'
import sys
p,n=sys.argv[1],sys.argv[2]; s=open(p).read()
old_bc='''! Bedrock
Boundary Condition 1
  Target Boundaries = 1
  Velocity 1 = Real 0.0e0
  Velocity 2 = Real 0.0e0
'''
new_bc='''! Bedrock
Boundary Condition 1
  Target Boundaries = 1
  Normal-Tangential Velocity = Logical True
  Flow Force BC = Logical True
  Velocity 1 = Real 0.0
  Slip Coefficient 2 =  Variable Coordinate 1
    Real Procedure "ElmerIceUSF" "%s"
'''%n
assert old_bc in s
s=s.replace(old_bc,new_bc,1)
s=s.replace('Solver 2\n  Equation = "Navier-Stokes"','Solver 2\n  Equation = "Normal vector"\n  Procedure = "ElmerIceSolvers" "ComputeNormalSolver"\n  Variable = -dofs 2 "Normal Vector"\n  Variable DOFs = 2\n  Optimize Bandwidth = Logical False\n  ComputeAll = Logical False\nEnd\n\nSolver 3\n  Equation = "Navier-Stokes"',1)
s=s.replace('Solver 3\n  Equation = Sij','Solver 4\n  Equation = Sij')
s=s.replace('Solver 4\n  Exec Solver = After All','Solver 5\n  Exec Solver = After All')
s=s.replace('Active Solvers(4) = 1 2 3 4','Active Solvers(5) = 1 2 3 4 5')
s=s.replace('Show Norm Index = 1','Show Norm Index = 1\n  Filename = scalars.dat')
open(p,'w').write(s)
PYEOF
}
mk_stress S-effectivepressure EffectivePressure
python3 - "$WORK/S-effectivepressure.sif" <<'PYEOF'
import sys
p=sys.argv[1]; s=open(p).read()
s=s.replace('Velocity 1 = Real 0.0\n','Velocity 1 = Real 0.0\n  External Pressure = Real 1.0e5\n',1)
open(p,'w').write(s)
PYEOF
run_pair S-effectivepressure
mk_stress S-gzept GZEPT_scaling
python3 - "$WORK/S-gzept.sif" <<'PYEOF'
import sys
p=sys.argv[1]; s=open(p).read()
s=s.replace('Velocity 1 = Real 0.0\n','Velocity 1 = Real 0.0\n  External Pressure = Real 1.0e5\n  GZEPT = Real 1.0e9\n',1)
open(p,'w').write(s)
PYEOF
run_pair S-gzept
mk_stress S-budd Sliding_Budd
python3 - "$WORK/S-budd.sif" <<'PYEOF'
import sys
p=sys.argv[1]; s=open(p).read()
s=s.replace('Velocity 1 = Real 0.0\n','Velocity 1 = Real 0.0\n  External Pressure = Real 1.0e5\n  Budd Friction Coefficient = Real 0.02\n  Budd Velocity Exponent = Real 3.0\n  Budd Zab Exponent = Real 1.0\n  Budd Linear Velocity = Real 1.0e-4\n  Budd Gravity = Real 9.81\n  Budd Floatation = Logical False\n',1)
open(p,'w').write(s)
PYEOF
run_pair S-budd
echo "STAGE-S3-GREEN $(date -u +%FT%TZ)"
