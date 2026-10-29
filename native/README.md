# native/

Rust-native Elmer rewrite directory, kept separate from the frozen ABI-port crates in `crates/`.

## Layout

```text
native/
  elmer-oracle/      # external oracle: direct runner, golden capture, grader (own workspace, own lockfile)
  docs/adr/          # stage decision records (oracle design, layout, arms)
  docker/            # pinned container definitions for the Fortran arms
  oracle-arms/       # per-arm configure/build scripts and version manifests
```

## Rules

- `elmer-oracle` is its own Cargo workspace with its own lockfile; a future `native/` swarm workspace explicitly excludes it, so candidate lanes can never become oracle dependencies.
- The oracle never depends on `crates/` or `tools/` (enforced by `elmer-oracle/tests/isolation.rs` on every `cargo test` run).
- New crates are licensed LGPL-2.1-or-later (decision D7); existing crates keep their licences.
- Golden data and build trees live outside the repo, under `/home/john/firstmate/data/elmer-s0-oracle/`.
