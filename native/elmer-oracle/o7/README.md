# O7 kernel fixtures (S0, report §5.2 item 8, §4 O7)

Oracle data. The oracle's Rust code calls the pristine Fortran kernels
directly and requires 0-ulp agreement. A whole-test mismatch then reduces
to a one-function diff.

| File | Content |
|---|---|
| `o7_shim.f90` | `BIND(C)` entry points over one arm's `libelmersolver`: `o7_rule` (the Integration.F90 rule builders), `o7_info` (`ElementInfo`: detJ, Basis, dBasisdx) and `o7_init`. |
| `rules.txt` | 89 Gauss rules: line 1–13, triangle, quad, tetra, pyramid, wedge, brick, `GaussPointsWedge2` and the economic wedge. `rule <kind> <args> <n>`, then `u v w s` per point. |
| `info.txt` | 792 `ElementInfo` calls on 51 curved or distorted elements of 21 Lagrange types (202 203 204 303 306 310 404 408 409 412 416 504 510 605 613 706 715 718 808 820 827) in 1D, 2D and 3D models. These are integration-point calls (`pt ip`) and the calls `StabParam` makes (`pt stab`). |

Reals are f64 bit patterns in hex, the format of the candidate element
library's captured data, so a candidate tests against these files without
depending on this crate. Every type S1 needs (§5.3: 202/303/404/504/808,
306/408/510/820, and 203 on their edges) is covered; `tests/o7.rs` checks
that on every `cargo test`.

## Provenance and result

The inputs (rule arguments, element types, nodes, call points) are those
of `elmer-fe/tests/data/{rules,info}.txt`, captured with gdb from the
pristine Debug build. `elmer-oracle o7-replay` re-evaluates every input on
each arm through the shim:

| Arm | Rule differences | `ElementInfo` differences |
|---|---:|---:|
| Release (golden) | **0** | **0** |
| OpenBLAS | 0 | 0 |
| FMA (`-march=x86-64-v3 -ffp-contract=fast`) | 3,287 | 13,900 |

The committed files are the Release outputs. They are byte-identical to
the Debug capture, so values the candidate matches at 0 ulp are also the
golden arm's. The FMA arm is the control: contraction changes the bits,
and the harness sees it.

## Reproduce

```sh
B=<arm build dir>            # e.g. .../elmer-s0-oracle/build-release
scripts/o7-build-shim.sh $B shim.so
target/release/elmer-oracle o7-replay $B shim.so o7 <out-dir>
ELMER_O7_BUILD=$B ELMER_O7_SHIM=$PWD/shim.so cargo test --release --test o7
```

`o7-replay` prints the difference counts and writes the arm's values to
`<out-dir>`. Without the two variables, the pristine-arm test in
`tests/o7.rs` reports itself skipped. The loader preloads the arm's own
libraries by absolute path, because `libelmersolver`'s RUNPATH names the
container build tree.
