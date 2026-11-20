//! Write the deterministic fixture matrices for the pristine-library harness:
//! `cargo run --example gen_fixtures -- <out-dir> [count]`.

#[path = "common/fixture_gen.rs"]
mod fixture_gen;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let dir = std::path::PathBuf::from(args.get(1).expect("usage: gen_fixtures <out-dir> [count]"));
    let count = args.get(2).map(|s| s.parse().expect("count")).unwrap_or(64);
    std::fs::create_dir_all(&dir).expect("mkdir");
    for f in fixture_gen::fixtures(count) {
        std::fs::write(
            dir.join(format!("{}.mat", f.name)),
            fixture_gen::to_text(&f),
        )
        .expect("write");
    }
}
