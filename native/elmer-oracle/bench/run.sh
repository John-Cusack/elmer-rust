#!/bin/bash
# S0 benchmark baselines (strategy report 4.5): Release Fortran Elmer with
# OpenMP, and with MPI, on scaled cases; median of the runs, with RSS.
#
#   run.sh --arm DIR --out DIR [--scale small|full] [--runs N]
#          [--mode omp|mpi] [--threads T] [--np P] [--cpus LIST] [--cases a,b]
#
# --arm   the oracle's Release build (oracle-builds.json `release`), used read-only
# --cpus  the cores to pin to (taskset / mpirun --cpu-set), e.g. 0-7
# Writes <out>/runs.tsv (one row per run) and <out>/summary.tsv (medians).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
arm= out= scale=small runs=5 mode=omp threads=8 np=4 cpus= cases=heat3d,elast3d,poisson2d,ns2d,eigen3d
while [ $# -gt 0 ]; do
  case $1 in
    --arm) arm=$2; shift 2 ;; --out) out=$2; shift 2 ;; --scale) scale=$2; shift 2 ;;
    --runs) runs=$2; shift 2 ;; --mode) mode=$2; shift 2 ;; --threads) threads=$2; shift 2 ;;
    --np) np=$2; shift 2 ;; --cpus) cpus=$2; shift 2 ;; --cases) cases=$2; shift 2 ;;
    *) echo "unknown argument $1" >&2; exit 2 ;;
  esac
done
[ -n "$arm" ] && [ -n "$out" ] || { echo "--arm and --out are required" >&2; exit 2; }
[ -n "$cpus" ] || cpus=0-$(( (mode == "omp" ? threads : np) - 1 ))

# Mesh divisions per case: full sizes from report 4.5, small for smoke tests.
size() {
  case "$scale:$1" in
    full:heat3d) echo 100 ;; full:elast3d) echo 60 ;; full:poisson2d) echo 2000 ;;
    full:ns2d) echo 300 ;; full:eigen3d) echo 30 ;;
    small:heat3d) echo 12 ;; small:elast3d) echo 8 ;; small:poisson2d) echo 60 ;;
    small:ns2d) echo 20 ;; small:eigen3d) echo 6 ;;
  esac
}
grid() { case $1 in poisson2d|ns2d) echo square2d ;; *) echo box3d ;; esac; }
dofs_per_node() { case $1 in elast3d|eigen3d|ns2d) echo 3 ;; *) echo 1 ;; esac; }
# ElmerGrid -partition for MPI: px py pz.
parts() { case $1 in poisson2d|ns2d) echo "2 $(( np / 2 )) 1" ;; *) echo "2 $(( np / 2 )) 1" ;; esac; }

# hwloc's GL component probes the X server and hangs without one.
export HWLOC_COMPONENTS=-gl
export ELMER_HOME=$arm/fem/src ELMER_LIB=$arm/fem/src/modules
export LD_LIBRARY_PATH=$arm/fem/src:$arm/matc/src:$arm/fhutiter/src:$arm/mathlibs/src/arpack:$arm/mathlibs/src/parpack${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=$([ "$mode" = omp ] && echo "$threads" || echo 1) OMP_PROC_BIND=true OMP_PLACES=cores
mkdir -p "$out"
[ -f "$out/runs.tsv" ] || printf 'case\tmode\tscale\tworkers\trun\tdofs\twall_s\telmer_real_s\tmax_rss_kb\n' > "$out/runs.tsv"

for c in ${cases//,/ }; do
  n=$(size "$c"); w="$out/$scale-$mode-$c"
  rm -rf "$w"; mkdir -p "$w"; cp "$here/cases/$c/case.sif" "$w/"
  sed "s/@N@/$n/g" "$here/cases/$(grid "$c").grd.in" > "$w/mesh.grd"
  (cd "$w" && "$arm/elmergrid/src/ElmerGrid" 1 2 mesh.grd > grid.log 2>&1)
  nodes=$(awk 'NR==1{print $1}' "$w/mesh/mesh.header")
  dofs=$(( nodes * $(dofs_per_node "$c") ))
  if [ "$mode" = mpi ]; then
    (cd "$w" && "$arm/elmergrid/src/ElmerGrid" 2 2 mesh -partition $(parts "$c") > part.log 2>&1)
  fi
  for r in $(seq 1 "$runs"); do
    rm -f "$w"/time.*
    if [ "$mode" = omp ]; then
      (cd "$w" && taskset -c "$cpus" /usr/bin/time -v -o time.0 "$arm/fem/src/ElmerSolver" case.sif > run-$r.log 2>&1)
      rss=$(awk -F': ' '/Maximum resident set size/{print $2}' "$w/time.0")
      wall=$(awk -F': ' '/Elapsed \(wall clock\)/{print $2}' "$w/time.0")
      workers=$threads
    else
      t0=$(date +%s.%N)
      (cd "$w" && mpirun --cpu-set "$cpus" --bind-to core -np "$np" \
        sh -c '/usr/bin/time -v -o time.$OMPI_COMM_WORLD_RANK "$0" case.sif' "$arm/fem/src/ElmerSolver_mpi" > run-$r.log 2>&1)
      t1=$(date +%s.%N)
      # Ranks' RSS summed; wall from the launcher.
      rss=$(cat "$w"/time.* | awk -F': ' '/Maximum resident set size/{s+=$2} END{print s}')
      wall=$(echo "$t1 - $t0" | bc -l)
      workers=$np
    fi
    wall=$(echo "$wall" | awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else if (NF==2) print $1*60+$2; else print $1 }')
    real=$(awk '/SOLVER TOTAL TIME\(CPU,REAL\)/{print $NF}' "$w/run-$r.log" | tail -1)
    grep -q "ELMER SOLVER FINISHED" "$w/run-$r.log" || { echo "$c run $r did not finish; see $w/run-$r.log" >&2; exit 1; }
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$c" "$mode" "$scale" "$workers" "$r" "$dofs" "$wall" "$real" "$rss" >> "$out/runs.tsv"
  done
done

# Medians per case, mode, scale and worker count.
python3 - "$out/runs.tsv" "$out/summary.tsv" <<'PY'
import csv, statistics, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
groups = {}
for r in rows:
    groups.setdefault((r["case"], r["mode"], r["scale"], r["workers"], r["dofs"]), []).append(r)
with open(sys.argv[2], "w") as f:
    f.write("case\tmode\tscale\tworkers\tdofs\truns\tmedian_wall_s\tmedian_elmer_real_s\tmedian_max_rss_kb\n")
    for k, rs in sorted(groups.items()):
        med = lambda key: statistics.median(float(r[key]) for r in rs if r[key])
        f.write("\t".join(k) + f"\t{len(rs)}\t{med('wall_s'):.3f}\t{med('elmer_real_s'):.3f}\t{med('max_rss_kb'):.0f}\n")
PY
column -t -s $'\t' "$out/summary.tsv"
