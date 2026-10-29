//! Reader for Elmer binary `.result` files (`Output File` with
//! `Binary Output = True`, written by `SaveResult` in `ModelDescription.F90`).
//!
//! Layout: ASCII header line `BINARY 3.<E>`, one NUL byte, then per saved
//! step a `Time:` block followed by one block per variable. Each variable
//! block is a NUL-terminated name, a Perm block (int32 size, int32 nonzero
//! count, then that many `(index, value)` int32 pairs; size 0 means NULL,
//! -1 with an int64 offset means same-as-previous), and the values as
//! native-endian doubles in Perm order. The host endianness marker is
//! checked; cross-endian files are refused rather than misread.

use crate::{io_with, OracleError, Result};
use std::path::Path;

/// One saved timestep.
#[derive(Debug, Clone, PartialEq)]
pub struct TimeStep {
    pub saves_done: i32,
    pub ntime: i32,
    pub sim_time: f64,
    pub vars: Vec<VariableData>,
}

/// One variable in one timestep. `values` are stored in Perm order, exactly
/// as written; node-space mapping needs the mesh and lives in the grader.
#[derive(Debug, Clone, PartialEq)]
pub struct VariableData {
    pub name: String,
    pub perm: Option<Vec<(i32, i32)>>,
    pub perm_size: i32,
    pub values: Vec<f64>,
}

/// A parsed binary `.result` file.
#[derive(Debug, Clone, PartialEq)]
pub struct ResultFile {
    pub steps: Vec<TimeStep>,
}

struct Cursor<'a> {
    buf: &'a [u8],
    pos: usize,
    path: String,
}

impl<'a> Cursor<'a> {
    fn eof(&self) -> bool {
        self.pos >= self.buf.len()
    }

    fn fail(&self, what: &str) -> OracleError {
        OracleError::Parse(format!("{}: {what} at byte {}", self.path, self.pos))
    }

    fn take(&mut self, n: usize, what: &str) -> Result<&'a [u8]> {
        if self.pos + n > self.buf.len() {
            return Err(self.fail(what));
        }
        let out = &self.buf[self.pos..self.pos + n];
        self.pos += n;
        Ok(out)
    }

    fn int4(&mut self) -> Result<i32> {
        let b = self.take(4, "int32")?;
        Ok(i32::from_ne_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn int8(&mut self) -> Result<i64> {
        let b = self.take(8, "int64")?;
        Ok(i64::from_ne_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ]))
    }

    fn double(&mut self) -> Result<f64> {
        let b = self.take(8, "double")?;
        Ok(f64::from_ne_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ]))
    }

    /// NUL-terminated string, as `BinWriteString` emits.
    fn cstring(&mut self) -> Result<String> {
        let start = self.pos;
        while self.pos < self.buf.len() && self.buf[self.pos] != 0 {
            self.pos += 1;
        }
        if self.pos >= self.buf.len() {
            return Err(self.fail("unterminated string"));
        }
        let s = std::str::from_utf8(&self.buf[start..self.pos])
            .map_err(|_| self.fail("non-UTF8 variable name"))?;
        self.pos += 1;
        Ok(s.to_string())
    }
}

fn host_endian_marker() -> u8 {
    if cfg!(target_endian = "little") {
        b'L'
    } else {
        b'B'
    }
}

/// Parse a binary `.result` file from memory.
pub fn parse_result_bytes(buf: &[u8], origin: &str) -> Result<ResultFile> {
    let mut cur = Cursor {
        buf,
        pos: 0,
        path: origin.to_string(),
    };
    // ASCII header: `BINARY 3.<E>` written list-directed, so leading blank.
    let header_end = buf
        .iter()
        .position(|&b| b == b'\n')
        .ok_or_else(|| cur.fail("missing header line"))?;
    let header = std::str::from_utf8(&buf[..header_end])
        .map_err(|_| cur.fail("non-UTF8 header"))?;
    let header = header.trim();
    if !header.starts_with("BINARY 3.") || header.len() != "BINARY 3.X".len() {
        return Err(cur.fail("not a binary result file"));
    }
    let marker = header.as_bytes()[header.len() - 1];
    if marker != host_endian_marker() {
        return Err(OracleError::Parse(format!(
            "{origin}: endian marker {marker} does not match host; refusing cross-endian read"
        )));
    }
    // `InitializeFile` writes an ASCII preamble (date, dof count, vector
    // section names) before `SwitchToBinary` appends the binary marker:
    // the first NUL byte ends the preamble.
    let nul_at = buf[header_end..]
        .iter()
        .position(|&b| b == 0)
        .ok_or_else(|| cur.fail("missing binary marker NUL"))?;
    cur.pos = header_end + nul_at + 1;

    let mut steps = Vec::new();
    let mut prev_offsets: Vec<(usize, usize)> = Vec::new();
    while !cur.eof() {
        let step_pos = cur.pos;
        let tag = cur.cstring()?;
        if tag != "Time:" {
            return Err(cur.fail(&format!("expected Time: block, found {tag:?}")));
        }
        let saves_done = cur.int4()?;
        let ntime = cur.int4()?;
        let sim_time = cur.double()?;
        let mut vars = Vec::new();
        // Variables run until the next `Time:` tag or EOF. Peek by scanning
        // for the next NUL-terminated string without consuming on failure:
        // instead, parse a var block, and stop when the upcoming tag parses
        // as `Time:` followed by a plausible step. Simpler and exact: the
        // count is unknown, so attempt var parse; on `Time:` tag, rewind.
        loop {
            if cur.eof() {
                break;
            }
            let probe = cur.pos;
            let name = cur.cstring()?;
            if name == "Time:" {
                cur.pos = probe;
                break;
            }
            let size = cur.int4()?;
            let (perm, perm_size): (Option<Vec<(i32, i32)>>, i32) = if size == 0 {
                (None, 0)
            } else if size == -1 {
                let off = cur.int8()?;
                let saved = cur.pos;
                if off < 0 || (off as usize) >= buf.len() {
                    return Err(cur.fail("same-as-previous offset out of range"));
                }
                cur.pos = off as usize;
                let psize = cur.int4()?;
                let pnnz = cur.int4()?;
                if psize < 0 || pnnz < 0 {
                    return Err(cur.fail("negative perm block at shared offset"));
                }
                let mut pairs = Vec::with_capacity(pnnz as usize);
                for _ in 0..pnnz {
                    pairs.push((cur.int4()?, cur.int4()?));
                }
                cur.pos = saved;
                prev_offsets.push((step_pos, pairs.len()));
                (Some(pairs), psize)
            } else {
                if size < 0 {
                    return Err(cur.fail("negative perm size"));
                }
                let nnz = cur.int4()?;
                if nnz < 0 || nnz > size {
                    return Err(cur.fail("perm count out of range"));
                }
                let mut pairs = Vec::with_capacity(nnz as usize);
                for _ in 0..nnz {
                    pairs.push((cur.int4()?, cur.int4()?));
                }
                (Some(pairs), size)
            };
            // Value count: variables with `DOFs == 1` store every Perm slot
            // with a positive entry; without a mesh the reader stores every
            // value that follows until the next name, which is ambiguous.
            // SaveResult writes exactly the positive Perm entries (or the
            // full array when Perm is NULL), so count from the block.
            let count = match &perm {
                None => {
                    return Err(cur.fail(
                        "NULL-Perm variables need mesh sizes; not supported yet",
                    ));
                }
                Some(pairs) => pairs.len(),
            };
            let mut values = Vec::with_capacity(count);
            for _ in 0..count {
                values.push(cur.double()?);
            }
            vars.push(VariableData {
                name,
                perm,
                perm_size,
                values,
            });
        }
        steps.push(TimeStep {
            saves_done,
            ntime,
            sim_time,
            vars,
        });
    }
    if steps.is_empty() {
        return Err(OracleError::Parse(format!("{origin}: no timesteps found")));
    }
    Ok(ResultFile { steps })
}

/// Read and parse a binary `.result` file from disk.
pub fn read_result_file(path: &Path) -> Result<ResultFile> {
    let buf = std::fs::read(path).map_err(|e| io_with(path, e))?;
    parse_result_bytes(&buf, &path.display().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn w_string(out: &mut Vec<u8>, s: &str) {
        out.extend_from_slice(s.as_bytes());
        out.push(0);
    }

    fn w_i32(out: &mut Vec<u8>, v: i32) {
        out.extend_from_slice(&v.to_ne_bytes());
    }

    fn w_i64(out: &mut Vec<u8>, v: i64) {
        out.extend_from_slice(&v.to_ne_bytes());
    }

    fn w_f64(out: &mut Vec<u8>, v: f64) {
        out.extend_from_slice(&v.to_ne_bytes());
    }

    fn sample_file_with(preamble: &[u8]) -> Vec<u8> {
        let mut b = Vec::new();
        b.extend_from_slice(b" BINARY 3.L\n");
        // Offset of the first Perm block: header + preamble + NUL +
        // Time tag + two int32s + one double + variable name.
        let perm_at = b.len() + preamble.len() + 1 + 6 + 4 + 4 + 8 + 12;
        b.extend_from_slice(preamble);
        b.push(0);
        w_string(&mut b, "Time:");
        w_i32(&mut b, 1);
        w_i32(&mut b, 1);
        w_f64(&mut b, 0.5);
        w_string(&mut b, "temperature");
        w_i32(&mut b, 3);
        w_i32(&mut b, 2);
        w_i32(&mut b, 1);
        w_i32(&mut b, 1);
        w_i32(&mut b, 3);
        w_i32(&mut b, 2);
        w_f64(&mut b, 300.0);
        w_f64(&mut b, 310.5);
        w_string(&mut b, "Time:");
        w_i32(&mut b, 2);
        w_i32(&mut b, 2);
        w_f64(&mut b, 1.0);
        w_string(&mut b, "temperature");
        w_i32(&mut b, -1);
        w_i64(&mut b, perm_at as i64);
        w_f64(&mut b, 301.0);
        w_f64(&mut b, 311.5);
        b
    }

    fn sample_file() -> Vec<u8> {
        sample_file_with(b"")
    }

    #[test]
    fn parses_two_steps_with_shared_perm() {
        let f = parse_result_bytes(&sample_file(), "test").unwrap();
        assert_eq!(f.steps.len(), 2);
        assert_eq!(f.steps[0].vars.len(), 1);
        let v0 = &f.steps[0].vars[0];
        assert_eq!(v0.name, "temperature");
        assert_eq!(v0.perm_size, 3);
        assert_eq!(v0.values, vec![300.0, 310.5]);
        // Second step reuses the first step's Perm via offset 24.
        let v1 = &f.steps[1].vars[0];
        assert_eq!(v1.perm, v0.perm);
        assert_eq!(v1.values, vec![301.0, 311.5]);
        assert_eq!(f.steps[1].sim_time, 1.0);
    }

    #[test]
    fn rejects_wrong_endian_marker() {
        let mut b = sample_file();
        let other = if cfg!(target_endian = "little") {
            b'B'
        } else {
            b'L'
        };
        let pos = b.iter().position(|&x| x == b'\n').unwrap();
        b[pos - 1] = other;
        assert!(parse_result_bytes(&b, "test").is_err());
    }

    #[test]
    fn rejects_truncated_file() {
        let b = sample_file();
        assert!(parse_result_bytes(&b[..b.len() - 3], "test").is_err());
    }

    #[test]
    fn rejects_garbage_header() {
        assert!(parse_result_bytes(b"hello\n\0", "test").is_err());
    }

    #[test]
    fn skips_ascii_preamble_before_binary_marker() {
        // Real files carry InitializeFile's text preamble between the
        // header line and the NUL marker.
        let file = sample_file_with(b"!File started at: 2026-10-29\n Degrees of freedom: \n");
        let f = parse_result_bytes(&file, "test").unwrap();
        assert_eq!(f.steps.len(), 2);
        assert_eq!(f.steps[0].vars[0].values, vec![300.0, 310.5]);
        assert_eq!(f.steps[1].vars[0].values, vec![301.0, 311.5]);
    }
}
