#!/bin/bash
# Grade a candidate ElmerSolver through the external grader (report §4.3
# item 1). Each listed test runs with the candidate's solver and the
# pristine arm's everything else (registration, ElmerGrid, staging), with
# the O2 injection, and every golden result file is then graded against
# the candidate's with the committed Tier bars. The candidate's own
# TEST.PASSED and PASS lines are ignored.
# Usage: grade-candidate.sh <elmer-oracle-bin> <candidate-solver> <build-dir> <tests-dir> <golden-root> <tolerances.json> <out> [jobs] < test-names
set -u
BIN="$1"; CAND="$2"; BUILD="$3"; TESTS="$4"; GOLDEN="$5"; TOL="$6"; OUT="$7"; JOBS="${8:-6}"
export HWLOC_COMPONENTS=-gl OMP_NUM_THREADS=1 OMPI_MCA_hwloc_base_binding_policy=none
export BIN CAND BUILD TESTS OUT
mkdir -p "$OUT"
# The candidate runs in a directory named like the test's own (as under
# ctest, `fem/tests/<dir>`): a candidate may key test-local fixtures by its
# working directory's name. `run` links to it, the layout grade-tree reads.
export GOLDEN
one() {
  dest="$OUT/$1"
  [ -e "$dest" ] && return
  mkdir -p "$dest"
  dir=$(awk -F'\t' -v t="$1" '$1 == t {print $2; exit}' "$GOLDEN/registry.tsv")
  dir=$(basename "${dir:-$1}")
  timeout 900 "$BIN" run "$BUILD" "$TESTS" "$1" "$dest/$dir" --o2 --solver "$CAND" > "$dest/run.log" 2>&1
  echo "$?" > "$dest/rc"
  [ "$dir" = run ] || ln -s "$dir" "$dest/run"
}
export -f one
xargs -P "$JOBS" -I{} bash -c 'one "$@"' _ {}
cp "$GOLDEN/registry.tsv" "$OUT/registry.tsv"
"$BIN" grade-tree "$TESTS" "$GOLDEN" "$OUT" "$TOL" > "$OUT/grades.tsv"
