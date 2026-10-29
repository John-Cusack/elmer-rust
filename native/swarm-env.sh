# Build environment for swarm lanes (swarm plan 4.1). Source it; do not run it:
#   . <worktree>/native/swarm-env.sh
# Builds: cargo <cmd> --manifest-path <worktree>/native/Cargo.toml -p <own-crate>
# Never: workspace-wide --release builds, cmake/make in elmer-fork, mpirun -np >1.
# Pristine runs: nice -n 10, serial ElmerSolver only, at most 2 concurrent per
# lane, in lane scratch.

export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-2}"   # L00 uses 4
export RUSTC_WRAPPER=sccache                      # shared cache across worktrees
export SCCACHE_CACHE_SIZE=30G
export CARGO_PROFILE_DEV_DEBUG=line-tables-only
export CARGO_PROFILE_TEST_DEBUG=line-tables-only
# No mold flag: the pinned nightly already links with rust-lld, and
# `-fuse-ld=mold` fails here (no ld.mold on PATH). Leave RUSTFLAGS unset so
# every lane shares one sccache key.
export HWLOC_COMPONENTS=-gl                       # any Elmer/MPI process
