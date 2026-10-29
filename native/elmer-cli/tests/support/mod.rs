//! Shared test helpers (tests only).
#![allow(dead_code)]

use std::io::Write;
use std::sync::{Arc, Mutex};

use elmer_cli::logger::{ElmerLogger, OutputConfig};

#[derive(Clone, Default)]
pub struct Buf(pub Arc<Mutex<Vec<u8>>>);

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
        String::from_utf8(self.0.lock().unwrap().clone()).unwrap()
    }
}

pub fn logger() -> (ElmerLogger, Buf) {
    let b = Buf::default();
    (ElmerLogger::new(OutputConfig::default(), Box::new(b.clone())), b)
}
