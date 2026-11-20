//! Rust side of the bitwise comparison harness: reads a matrix file and
//! prints the same lines as the pristine-library harness.
//!
//! `cargo run --release --example umf_harness -- amd <matrix>`

#[path = "common/harness_io.rs"]
mod harness_io;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert!(args.len() >= 3, "usage: umf_harness <mode> <matrix>");
    let text = std::fs::read_to_string(&args[2]).expect("read matrix");
    print!("{}", harness_io::run(&args[1], &harness_io::parse(&text)));
}
