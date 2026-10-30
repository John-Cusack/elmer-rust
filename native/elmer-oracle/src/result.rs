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
    if header == "ASCII 3" {
        return parse_ascii_result(buf, origin);
    }
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
    let sizes = preamble_sizes(&buf[header_end..header_end + nul_at]);
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
                None => null_perm_count(&sizes, &name).ok_or_else(|| {
                    cur.fail(&format!("NULL-Perm variable {name:?} missing from the preamble"))
                })?,
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

/// `SIZE(Var % Values)` per variable from the `Degrees of freedom:`
/// preamble `InitializeFile` writes (`name : size permsize dofs : eq`). A
/// variable with a NULL Perm stores its whole value array, and this is the
/// only place the file records how long that is.
fn preamble_sizes(preamble: &[u8]) -> Vec<(String, usize)> {
    let text = String::from_utf8_lossy(preamble);
    let mut out = Vec::new();
    for line in text.lines() {
        let mut parts = line.splitn(3, " : ");
        let (Some(name), Some(nums)) = (parts.next(), parts.next()) else {
            continue;
        };
        if let Some(size) = nums.split_whitespace().next().and_then(|n| n.parse().ok()) {
            out.push((name.trim().to_string(), size));
        }
    }
    out
}

fn null_perm_count(sizes: &[(String, usize)], name: &str) -> Option<usize> {
    sizes.iter().find(|(n, _)| n == name).map(|&(_, s)| s)
}

/// Parse the ASCII form (`Binary Output` off): the same blocks as text,
/// values list-directed (17 significant digits, so the round trip is
/// exact), `Perm: use previous` meaning the last Perm written.
fn parse_ascii_result(buf: &[u8], origin: &str) -> Result<ResultFile> {
    let fail = |what: &str| OracleError::Parse(format!("{origin}: {what}"));
    let text = String::from_utf8_lossy(buf);
    let mut lines = text.lines().peekable();
    let mut preamble = String::new();
    while let Some(l) = lines.peek() {
        if l.starts_with("Time:") {
            break;
        }
        preamble.push_str(l);
        preamble.push('\n');
        lines.next();
    }
    let sizes = preamble_sizes(preamble.as_bytes());
    let mut steps = Vec::new();
    let mut prev_perm: Option<(Vec<(i32, i32)>, i32)> = None;
    while let Some(line) = lines.next() {
        let rest = line
            .strip_prefix("Time:")
            .ok_or_else(|| fail(&format!("expected Time: line, found {line:?}")))?;
        let mut nums = rest.split_whitespace();
        let mut next_num = |what: &str| nums.next().ok_or_else(|| fail(what));
        let saves_done = next_num("Time: saves")?.parse().map_err(|_| fail("Time: saves"))?;
        let ntime = next_num("Time: step")?.parse().map_err(|_| fail("Time: step"))?;
        let sim_time = next_num("Time: time")?.parse().map_err(|_| fail("Time: time"))?;
        let mut vars = Vec::new();
        while let Some(l) = lines.peek() {
            if l.starts_with("Time:") {
                break;
            }
            let name = l.trim_end().to_string();
            lines.next();
            let perm_line = lines.next().ok_or_else(|| fail("missing Perm line"))?;
            let perm_spec = perm_line
                .strip_prefix("Perm:")
                .ok_or_else(|| fail(&format!("expected Perm: for {name:?}")))?
                .trim();
            let (perm, perm_size) = match perm_spec {
                "NULL" => (None, 0),
                "use previous" => {
                    let (pairs, size) =
                        prev_perm.clone().ok_or_else(|| fail("use previous with no Perm yet"))?;
                    (Some(pairs), size)
                }
                spec => {
                    let mut it = spec.split_whitespace().map(|n| n.parse::<i32>());
                    let size = it.next().and_then(|r| r.ok()).ok_or_else(|| fail("Perm size"))?;
                    let nnz = it.next().and_then(|r| r.ok()).ok_or_else(|| fail("Perm count"))?;
                    let mut pairs = Vec::with_capacity(nnz.max(0) as usize);
                    for _ in 0..nnz {
                        let pl = lines.next().ok_or_else(|| fail("truncated Perm"))?;
                        let mut p = pl.split_whitespace().map(|n| n.parse::<i32>());
                        match (p.next(), p.next()) {
                            (Some(Ok(a)), Some(Ok(b))) => pairs.push((a, b)),
                            _ => return Err(fail(&format!("bad Perm pair {pl:?}"))),
                        }
                    }
                    prev_perm = Some((pairs.clone(), size));
                    (Some(pairs), size)
                }
            };
            let count = match &perm {
                None => null_perm_count(&sizes, &name)
                    .ok_or_else(|| fail(&format!("NULL-Perm variable {name:?} not in preamble")))?,
                Some(pairs) => pairs.len(),
            };
            let mut values = Vec::with_capacity(count);
            for _ in 0..count {
                let vl = lines.next().ok_or_else(|| fail("truncated values"))?;
                values.push(
                    vl.trim()
                        .parse::<f64>()
                        .map_err(|_| fail(&format!("bad value {vl:?}")))?,
                );
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
        return Err(fail("no timesteps found"));
    }
    Ok(ResultFile { steps })
}

/// `SIZE(Var % Values)` per variable name, from a result file's preamble
/// (either form): the full length of each solver array, whether or not
/// every entry is reachable through the Perm.
pub fn value_sizes(buf: &[u8]) -> Vec<(String, usize)> {
    let end = buf
        .iter()
        .position(|&b| b == 0)
        .or_else(|| buf.windows(6).position(|w| w == b"\nTime:"))
        .unwrap_or(buf.len());
    preamble_sizes(&buf[..end])
}

/// `Number Of Nodes:` from a result file's preamble (the mesh's nodal
/// count; a longer Perm carries edge, face or DG entries).
pub fn number_of_nodes(buf: &[u8]) -> Option<usize> {
    let head = &buf[..buf.len().min(1 << 16)];
    let text = String::from_utf8_lossy(head);
    text.lines()
        .find_map(|l| l.trim().strip_prefix("Number Of Nodes:"))
        .and_then(|n| n.trim().parse().ok())
}

/// Read and parse a `.result` file (binary or ASCII) from disk.
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
    fn ascii_form_parses_with_use_previous_and_null_perm() {
        let text = " ASCII 3\n!File started at: x\n Degrees of freedom: \n\
temperature              :        3       3   1 : heat equation\n\
cost                     :        2       0   1 : heat equation\n\
 Total DOFs:            1\n\
Time:       1      1  1.00000000E+000\n\
temperature\nPerm:            3            2\n          1          2\n          3          1\n\
   1.5000000000000000     \n   2.2500000000000000     \n\
cost\nPerm: NULL\n   7.0000000000000000     \n   8.0000000000000000     \n\
Time:       2      2  2.00000000E+000\n\
temperature\nPerm: use previous\n   3.0000000000000000     \n   4.0000000000000000     \n";
        let r = parse_result_bytes(text.as_bytes(), "t").unwrap();
        assert_eq!(r.steps.len(), 2);
        assert_eq!(r.steps[0].vars[0].perm, Some(vec![(1, 2), (3, 1)]));
        assert_eq!(r.steps[0].vars[0].values, vec![1.5, 2.25]);
        assert_eq!(r.steps[0].vars[1].perm, None);
        assert_eq!(r.steps[0].vars[1].values, vec![7.0, 8.0]);
        assert_eq!(r.steps[1].vars[0].perm, Some(vec![(1, 2), (3, 1)]));
        assert_eq!(r.steps[1].vars[0].values, vec![3.0, 4.0]);
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
