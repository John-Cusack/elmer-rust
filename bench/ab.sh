#!/bin/sh
# ab.sh: run the A/B measurement for one ported unit.
#
# usage:
#   bench/ab.sh <before-tree> <after-tree> <crate-dir> <bench-target> [bench args...]
#
#   before-tree   checkout of the port BEFORE the recraft (same bench files)
#   after-tree    checkout of the port AFTER the recraft
#   crate-dir     path of the port crate, relative to each tree
#                 (e.g. crates/elmer-fem-a/ports/ludecomp)
#   bench-target  the `[[bench]]` target name (e.g. ludecomp_bench)
#   bench args    passed through to the bench binary (e.g. --n 128 --iters 5)
#
# The script builds both arms, then lets the shared harness drive the A/B run
# (rep-interleaved, fresh child per sample, median + MAD, wait4 peak RSS,
# digest equality). Results land in <after-tree>/<crate-dir>/bench-out/.
#
# Nothing here decides a verdict: the harness does, against its noise band.
set -eu

if [ "$#" -lt 4 ]; then
	echo "usage: $0 <before-tree> <after-tree> <crate-dir> <bench-target> [args...]" >&2
	exit 2
fi

before_tree="$1"; after_tree="$2"; crate_dir="$3"; bench_target="$4"; shift 4

# Build one arm and print the bench binary path.
build_arm() {
	tree="$1"
	(cd "$tree/$crate_dir" && cargo build --release --bench "$bench_target" >&2)
	# cargo puts `harness = false` bench binaries under
	# <workspace-root>/target/release/build/<pkg>/<hash>/out/<target>-<hash>.
	# The workspace root is not necessarily the crate dir (nested lane
	# workspaces), so walk up until the target dir shows up.
	dir="$tree/$crate_dir"
	while [ "$dir" != "/" ]; do
		if [ -d "$dir/target/release/build" ]; then
			bin=$(ls -t "$dir/target/release/build"/*/*/out/"$bench_target"-* 2>/dev/null \
				| grep -v '\.d$' | head -1)
			if [ -n "$bin" ]; then
				echo "$bin"
				return 0
			fi
		fi
		dir=$(dirname "$dir")
	done
	return 1
}

before_bin="$(build_arm "$before_tree")"
after_bin="$(build_arm "$after_tree")"
if [ -z "$before_bin" ] || [ -z "$after_bin" ]; then
	echo "ab.sh: could not locate a bench binary for '$bench_target'" >&2
	exit 1
fi
echo "ab.sh: before=$before_bin" >&2
echo "ab.sh: after=$after_bin" >&2

out_dir=""
prev=""
for arg in "$@"; do
	if [ "$prev" = "--out" ]; then
		out_dir="$arg"
	fi
	prev="$arg"
done
if [ -z "$out_dir" ]; then
	out_dir="$after_tree/bench-out/$bench_target"
fi
mkdir -p "$out_dir"
# shellcheck disable=SC2086
"$after_bin" --mode parent \
	--before-bin "$before_bin" \
	--after-bin "$after_bin" \
	--out "$out_dir" \
	"$@"
