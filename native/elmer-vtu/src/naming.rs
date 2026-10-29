//! File names, `VtuFileNaming` (VtuOutputSolver.F90:104-190).

/// Everything `VtuFileNaming` reads besides its arguments.
#[derive(Clone, Debug, PartialEq)]
pub struct VtuNaming {
    /// `ParEnv % PEs`.
    pub pes: usize,
    /// `No Fileindex` of the solver section.
    pub no_fileindex: bool,
    /// `Partition Numbering`: tag wrapper files (no `part`) with `<PEs>np`.
    pub parallel_base: bool,
    /// Group name (`_Body1`, `_BC2`, a section `Name`) for `Vtu Part Collection`.
    pub group_name: Option<String>,
}

impl Default for VtuNaming {
    fn default() -> Self {
        VtuNaming { pes: 1, no_fileindex: false, parallel_base: false, group_name: None }
    }
}

/// Zero-padded width for the partition count in a name (`I1.1` .. `I4.4`).
fn pe_width(pes: usize) -> usize {
    match pes {
        0..=9 => 1,
        10..=99 => 2,
        100..=999 => 3,
        _ => 4,
    }
}

impl VtuNaming {
    /// Compose `base` + group + partition + time index + `suffix`.
    ///
    /// * `part` is the 1-based partition of a `.vtu` piece; wrapper files
    ///   (`.pvtu`, `.pvd`) pass `None`.
    /// * `file_index <= 0` or `no_fileindex` omits the `_tNNNN` part, which is
    ///   four digits below 10000 and as many as needed above.
    /// * `no_path` drops everything up to the last `/` (names stored inside
    ///   `.pvd` / `.pvtu` are relative to the wrapper).
    ///
    /// With one partition the partition tag is never written.
    pub fn file_name(
        &self,
        base: &str,
        suffix: &str,
        file_index: i32,
        part: Option<usize>,
        no_path: bool,
    ) -> String {
        let mut name = base.to_string();

        if let Some(g) = &self.group_name {
            name.push('_');
            name.push_str(g);
        }

        if self.pes != 1 {
            let w = pe_width(self.pes);
            match part {
                Some(p) => name.push_str(&format!("_{:0w$}np{:0w$}", self.pes, p, w = w)),
                None if self.parallel_base => name.push_str(&format!("_{:0w$}np", self.pes, w = w)),
                None => {}
            }
        }

        if file_index > 0 && !self.no_fileindex {
            if file_index < 10000 {
                name.push_str(&format!("_t{file_index:04}"));
            } else {
                name.push_str(&format!("_t{file_index}"));
            }
        }

        if no_path {
            if let Some(j) = name.rfind('/') {
                name.drain(..=j);
            }
        }

        name.push_str(suffix);
        name
    }
}

/// `case_t0001.vtu`-style name for a run with `pes` partitions (no groups,
/// no `Partition Numbering`): `base[_<pes>np<part>][_tNNNN]<suffix>`.
pub fn vtu_file_name(
    base: &str,
    suffix: &str,
    file_index: i32,
    part: Option<usize>,
    pes: usize,
) -> String {
    VtuNaming { pes, ..VtuNaming::default() }.file_name(base, suffix, file_index, part, false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn serial_names() {
        assert_eq!(vtu_file_name("out/case", ".vtu", 1, Some(1), 1), "out/case_t0001.vtu");
        assert_eq!(vtu_file_name("case", ".vtu", 12345, Some(1), 1), "case_t12345.vtu");
        assert_eq!(vtu_file_name("case", ".vtu", 0, Some(1), 1), "case.vtu");
    }

    #[test]
    fn partition_widths() {
        assert_eq!(vtu_file_name("c", ".vtu", 1, Some(3), 4), "c_4np3_t0001.vtu");
        assert_eq!(vtu_file_name("c", ".vtu", 1, Some(3), 12), "c_12np03_t0001.vtu");
        assert_eq!(vtu_file_name("c", ".vtu", 1, Some(7), 128), "c_128np007_t0001.vtu");
        assert_eq!(vtu_file_name("c", ".vtu", 1, Some(7), 1024), "c_1024np0007_t0001.vtu");
        // Wrappers carry no partition unless Partition Numbering is on.
        assert_eq!(vtu_file_name("c", ".pvtu", 1, None, 4), "c_t0001.pvtu");
        let n = VtuNaming { pes: 4, parallel_base: true, ..Default::default() };
        assert_eq!(n.file_name("c", ".pvtu", 1, None, false), "c_4np_t0001.pvtu");
    }

    #[test]
    fn groups_no_fileindex_and_no_path() {
        let n = VtuNaming {
            pes: 2,
            group_name: Some("Body1".into()),
            no_fileindex: true,
            ..Default::default()
        };
        assert_eq!(n.file_name("d/c", ".vtu", 5, Some(2), true), "c_Body1_2np2.vtu");
    }
}
