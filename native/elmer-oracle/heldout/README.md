# Held-out seal (S0, report §4.1 O5 and §5.2 item 9)

Oracle data. Resealing needs an ADR with captain sign-off; the script
refuses to overwrite an existing seal.

## What is committed and what is not

- `seal.json` (committed): the method, the pool, the per-stratum counts, the
  exclusions, and `commitment_sha256`, the sha256 of the secret
  `sealed.json`.
- **`sealed.json` (not committed).** It holds the sealed directory names and
  the salt. It lives only with the oracle owner, at
  `/home/john/firstmate/data/elmer-s0-oracle-opus/heldout/secret/sealed.json`
  (mode 600). Anyone can check that file against the commitment:
  `sha256sum sealed.json` must print
  `7ddcaa627e895d18fbc641133bc3d8f0f3c461d50ba1541c32fc0e1a626cf4dd`.
- **The lane tree** is `fem/tests` and `elmerice/Tests` at the pin
  `9f6af2f85`, with the sealed directories removed. It is at
  `/home/john/firstmate/data/elmer-s0-oracle-opus/heldout/lane-tree`. Its
  manifest (`elmer-oracle manifest`, 6,694 files) has sha256
  `e6e1dcfc8a8baba667b985b3618830e8d64c677f17c7d811cfcbd9b2ec7a3c0a`.
  Candidate lanes get this tree, never the full one.

## Method

- **Unit: a test directory.** A directory that any quick test uses is never
  sealed, because the lane tree removes whole directories and the golden
  corpus needs every quick directory.
- **Pool: 520 directories.** These are the 473 `fem/tests` directories with
  no quick test, plus 47 `elmerice/Tests` directories. Excluded:
  - 76 directories that compile a test-local Fortran plugin. Decision D12
    keeps them out of held-out, and that includes 15 of the 65 ice tests.
  - 7 directories with a test that failed or timed out in the pristine
    ctest baseline (11 failures out of 683 non-quick tests). A test that
    fails in pristine cannot be graded against a pass.
- **Strata.** fem directories are stratified by serial or parallel, crossed
  with their most common family label (`whitney`, `mortar`, `mgdyn`,
  `shell`, ...), which gives 47 strata. 30% of each stratum is sealed,
  rounded half up. Directories are ranked by `sha256(salt + name)`, where the
  salt is a secret 128-bit random value.
- **Ice.** Every eligible ice directory is sealed, as §5.2 asks.

**Sealed: 190 directories, covering 204 registered tests.** That is 143 fem
directories plus all 47 ice directories, above the minimum of 21 (a
10-point gap must mean more than about 2 tests).

## Limits

- **Single-directory strata are never sealed** (11 of the 47 fem strata). Their
  families have no held-out signal.
- **The upstream source at the pin still contains every test.** The seal
  holds only as long as lanes are graded from, and pointed at, the lane
  tree, and do not browse upstream `fem/tests`. This is a process rule for
  firstmate, not something the oracle can enforce.
- **Item 9's first bullet is half done.** The pass/fail half of the
  non-quick baseline comes from the pristine ctest run. Pristine-vs-pristine
  determinism, which needs two captured runs of all 618 non-quick tests
  (about 70 min each at 8 jobs), has not been run yet.

## Reproduce

```sh
export ELMER_ORACLE=target/release/elmer-oracle
S=/home/john/firstmate/data/elmer-s0-oracle
scripts/heldout-seal.py seal $S/build-release $S/src $S/nonquick-baseline.log <secret-dir> heldout
scripts/heldout-seal.py lane-tree $S/src <secret-dir> <out>
```

A rerun draws a new salt and therefore a different seal. That is why
`seal` refuses to overwrite an existing one.
