//! UMFPACK numeric factorisation and solve, bitwise against the pristine
//! Elmer-bundled UMFPACK 4.4 on Elmer's path (`umfpack_di_defaults`,
//! `umfpack_di_symbolic`, `umfpack_di_numeric`, `umfpack_di_solve` with
//! `sys = UMFPACK_Aat` and no iterative refinement, as `umf4sol` calls it):
//! the symbolic lines, every numeric and solve Info entry (timings and the
//! allocator-dependent `UMFPACK_NUMERIC_COSTLY_REALLOC` excepted) and the
//! solution bits. Expected outputs (line hashes) come from the C library.

#[path = "../examples/common/fixture_gen.rs"]
mod fixture_gen;
#[path = "../examples/common/harness_io.rs"]
mod harness_io;

mod support;

#[test]
fn numeric_and_solve_match_pristine_bitwise() {
    let expected = support::expected("numeric_expected.txt");
    let mut checked = 0;
    for f in fixture_gen::fixtures(64) {
        let m = harness_io::parse(&fixture_gen::to_text(&f));
        support::check_hashed(&expected, &f.name, &harness_io::run("sol", &m));
        checked += 1;
    }
    for name in ["bentonite", "ElasticAnisotropic2D", "ElasticStressLoad"] {
        let m = harness_io::parse(&support::data(&format!("elmer_{name}.mat")));
        support::check_hashed(
            &expected,
            &format!("elmer_{name}"),
            &harness_io::run("sol", &m),
        );
        checked += 1;
    }
    assert_eq!(checked, expected.len());
}
