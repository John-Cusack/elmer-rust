//! Test-only stand-ins for S1's `elmer-mesh` (mesh-file reader), `elmer-fe`
//! (basis, quadrature, copied from `elmer-physics-heat`'s tests) and
//! `elmer-sif` (the hand-built MinimalExample model). Replace them with the
//! real crates once S1 merges.
#![allow(dead_code, clippy::needless_range_loop)]

pub mod basis;
pub mod case;
pub mod mesh;
pub mod quad;

use std::io::Write;
use std::sync::{Arc, Mutex};

use elmer_cli::logger::OutputConfig;
use elmer_cli::ElmerLogger;

/// A logger whose output the test can read back.
#[derive(Clone, Default)]
pub struct Buf(Arc<Mutex<Vec<u8>>>);

impl Write for Buf {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl Buf {
    pub fn text(&self) -> String {
        String::from_utf8_lossy(&self.0.lock().unwrap()).into_owned()
    }
}

pub fn logger() -> (ElmerLogger, Buf) {
    let buf = Buf::default();
    (ElmerLogger::new(OutputConfig::default(), Box::new(buf.clone())), buf)
}
