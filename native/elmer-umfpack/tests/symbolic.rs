//! UMFPACK symbolic analysis (`umfpack_di_symbolic`, default Control), bitwise
//! against the pristine Elmer-bundled UMFPACK 4.4: every Info entry (timings
//! excepted), the `umfpack_get_symbolic` arrays and the internal Symbolic
//! fields the numeric factorisation reads. Expected outputs (line hashes)
//! come from the C library on the same matrices; the fixtures cover the
//! symmetric, 2-by-2 and unsymmetric (COLAMD) strategies.

#[path = "../examples/common/fixture_gen.rs"]
mod fixture_gen;
#[path = "../examples/common/harness_io.rs"]
mod harness_io;

mod support;

#[test]
fn symbolic_matches_pristine_bitwise() {
    let expected = support::expected("symbolic_expected.txt");
    let mut checked = 0;
    for f in fixture_gen::fixtures(64) {
        let m = harness_io::parse(&fixture_gen::to_text(&f));
        support::check_hashed(&expected, &f.name, &harness_io::run("sym", &m));
        checked += 1;
    }
    for name in ["bentonite", "ElasticAnisotropic2D", "ElasticStressLoad"] {
        let m = harness_io::parse(&support::data(&format!("elmer_{name}.mat")));
        support::check_hashed(
            &expected,
            &format!("elmer_{name}"),
            &harness_io::run("sym", &m),
        );
        checked += 1;
    }
    assert_eq!(checked, expected.len());
}
