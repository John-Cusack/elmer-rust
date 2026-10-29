//! `Poisson` `PoissonSolver` with the 1dtests `Source` end to end against
//! the pristine ElmerSolver (26.2-devel, Rev 9f6af2f85) runs of
//! `fem/tests/1dtests` and `1dtests_MeshLevels` with the compiled
//! `Poisson.so` and `Source.so` (`ElmerGrid 1 2 1d -decimals 20`). Five
//! transient steps; `Source` checks each step's solution against its
//! analytic form and stops the run if it is off.

mod case;

const CASES: [(&str, &str); 2] = [("1dtests", "1.40257363E+00"), ("1dtests_MeshLevels", "1.41393619E+00")];

#[test]
fn one_d_tests_match_pristine() {
    for (name, pristine) in CASES {
        let dir = case::stage(name, None);
        let out = case::run(name, &dir, "1d.sif").unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(out.passed, Some(true), "{name}: reference check\n{}", out.log);
        case::assert_norm_prints(name, out.norm(1), pristine);
    }
}
