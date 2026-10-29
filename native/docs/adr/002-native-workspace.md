# ADR-002: layout of the Rust-native rewrite (`native/`)

Date: 2026-10-29. Slice S1 front-end (`elmer-sif`, `elmer-mesh`, `elmer-fe`).
(ADR-001 is owned by the oracle lane.)
workspace ("new repo", §10 D3). The captain chose the existing
John-Cusack/elmer-rust repo over a new one, with the native rewrite in its
own Cargo workspace directory kept separate from the frozen ABI-port crates,
without modifying or depending on them.

## Decision

- `native/` is a **self-contained Cargo workspace**. The workspace
  manifest `native/Cargo.toml` is owned by the skeleton lane
  (`fm/elmer-swarm-l00`, `members = ["elmer-*"]`); this lane does not
  commit it. A local copy is kept for builds until the rebase.
- No native crate depends on any `crates/*`, `fem-modules/*` or
  `tools/*` crate. Cross-lane integration goes through the shared
  `elmer-core` crate (skeleton lane): `elmer-mesh` implements
  `elmer_core::MeshView`, `elmer-fe` implements
  `elmer_core::{Quadrature, ElementBasis}`, and `elmer-sif` converts a
  section losslessly into `elmer_core::ParamMap`/`Value`. Thin adapters
  over this lane's own types; internal design unchanged.
- The top-level `Cargo.toml` is **not modified** by this lane.
- Shared toolchain: `../rust-toolchain.toml` still applies (no nested
  toolchain file), so the whole repo builds with one compiler.
- Snapshots of Elmer inputs consumed at build time live next to the
  crate that consumes them and record their origin pin:
  - `native/elmer-sif/keywords/SOLVER.KEYWORDS` (upstream `9f6af2f85`)
  - `native/elmer-mesh/elements.def` (upstream `9f6af2f85`)
  - `native/elmer-sif/matc/` (C sources from `matc/src`, LGPL-2.1)
- New crates are `LGPL-2.1-or-later`, each with its own `LICENSE` file;
  existing crates' licences are untouched (decision D7).

## Consequences

- `cargo build` at the repo root never touches `native/`; agents working
  in `native/` run cargo with `native/` as CWD (local manifest until
  the skeleton lands).
- Zero-Fortran enforcement (`native/check-zero-fortran.sh`, owned by
  this lane) scans only the `native/` workspace lockfile/tree, so frozen
  ABI-port crates (which legitimately invoke gfortran) cannot trip it.
- MATC C sources are compiled from the vendored snapshot by `cc`, so the
  `elmer-sif` build needs no Elmer tree and no Fortran compiler.
