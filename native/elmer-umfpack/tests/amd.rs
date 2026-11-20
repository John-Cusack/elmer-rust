//! AMD ordering, bitwise against the pristine Elmer-bundled AMD 1.1
//! (`amd_order`, default Control): permutation and every Info entry.
//! Expected outputs were produced by the C library from build-release on the
//! same matrices (the generated fixtures and three systems Elmer saved with
//! `Linear System Save` in ElasticAnisotropic2D, ElasticStressLoad, bentonite).

#[path = "../examples/common/fixture_gen.rs"]
mod fixture_gen;
#[path = "../examples/common/harness_io.rs"]
mod harness_io;

mod support;

#[test]
fn amd_order_matches_pristine_bitwise() {
    let expected = support::expected("amd_expected.txt");
    let mut checked = 0;
    for f in fixture_gen::fixtures(64) {
        let m = harness_io::parse(&fixture_gen::to_text(&f));
        support::check(&expected, &f.name, &harness_io::run("amd", &m));
        checked += 1;
    }
    for name in ["bentonite", "ElasticAnisotropic2D", "ElasticStressLoad"] {
        let m = harness_io::parse(&support::data(&format!("elmer_{name}.mat")));
        support::check(
            &expected,
            &format!("elmer_{name}"),
            &harness_io::run("amd", &m),
        );
        checked += 1;
    }
    assert_eq!(checked, expected.len());
}
