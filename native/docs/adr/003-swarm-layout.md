# ADR 003: swarm workspace layout

- Status: accepted.
- Date: 2026-10-29.
- Decides: the `native/` Cargo workspace shape, oracle isolation from candidate crates, and directory ownership for the parallel crate lanes (swarm plan sections 2.1, 2.2 and 4).

## Context

- The Rust-native rewrite is built by up to 16 concurrent lanes, each producing one crate under `native/`.
- If every lane edits the `members` line of a shared `native/Cargo.toml`, every merge after the first one conflicts.
- The oracle (`elmer-oracle`, ADR 001) must not share a dependency resolution with any crate it grades. Its isolation test rejects any candidate `elmer-*` package in its lockfile.
- Lanes compile against each other's types before any of them merge, so the types they share need a single owner that merges first.

## Decision: glob members

- `native/Cargo.toml` declares `members = ["elmer-*"]`. Adding a crate directory makes it a member; no lane edits the workspace manifest.
- `[workspace.package]` holds `edition = "2021"`, `license = "LGPL-2.1-or-later"` and the repository URL. Member crates inherit them.
- `[workspace.dependencies]` is the pre-resolved external allowlist, pinned in `native/Cargo.lock`. `elmer-core` pins the full allowlist through a `cfg(any())` dev-dependency table: cargo locks those entries but never builds them. A new external crate must be named in the PR description and needs firstmate's OK.
- `native/Cargo.lock` conflicts between lanes are mechanical. To resolve one, check out `main`'s lockfile and rerun `cargo metadata`, which adds the new path crates and keeps the existing pins. Never hand-edit it, and never run `cargo update` without `-p`.

## Decision: oracle as a separate workspace

- `native/Cargo.toml` sets `exclude = ["elmer-oracle"]`.
- `native/elmer-oracle/` is its own workspace root with its own `Cargo.lock`, and the oracle lane builds it with `--manifest-path native/elmer-oracle/Cargo.toml`.
- Rationale: the grader's lockfile can then only contain the grader, which is a stronger guarantee than a shared lockfile with a filter.

## Decision: one owner per directory

- Each lane owns exactly one `native/<crate>/` directory and writes nowhere else in `native/`.
- `native/Cargo.toml`, the initial `native/Cargo.lock`, `native/elmer-core/`, `native/docs/interfaces.md`, this ADR and `native/swarm-env.sh` belong to the skeleton lane (L00). After L00 merges, firstmate gatekeeps `elmer-core` through `core-change` PRs (`native/docs/interfaces.md`).
- `native/elmer-{sif,mesh,fe}/` and `native/check-zero-fortran.sh` belong to the front-end lane (S1). `native/elmer-oracle/` belongs to the oracle lane (S0).
- ADR numbers are assigned by firstmate: 001 oracle, 002 front-end, 003 this layout.

## Decision: `elmer-core` is small and frozen

- `elmer-core` holds only types that cross a lane boundary: the scalar trait, the CSR struct, the linear-operator, preconditioner and solver traits, the communicator, the parameter view, the mesh, basis, field and local-system views, and the logger.
- It has no logic beyond trivial constructors, `SerialComm`, `ParamMap`, `IdentityPrecond`, and the serial CSR product that only this crate is allowed to implement.
- Each lane keeps its own types in its own crate. `native/docs/interfaces.md` records who owns, implements and consumes each item.

## Build environment

- Lanes source `native/swarm-env.sh` and always run cargo with `--manifest-path <worktree>/native/Cargo.toml -p <crate>`. Firstmate runs the only workspace-wide build, serially, after each merge.
- The plan's `-fuse-ld=mold` link flag is omitted: no `ld.mold` is installed on this host (so `cc -fuse-ld=mold` fails), and the pinned nightly already links with rust-lld.
