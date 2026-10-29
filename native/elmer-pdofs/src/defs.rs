//! Element definitions: the `Element = "..."` keyword and the
//! `Solver % Def_Dofs(10, bodies, 6)` table it fills.
//!
//! Ports `GetDefs` (MeshBasics.F90:76), which `LoadModel`
//! (ModelDescription.F90:2883) runs on the solver-section keyword, and
//! `GetMaxDefs` (MeshBasics.F90:187), which `NonNodalElements`
//! (MeshLoad.F90:713) runs for body-wise and MATC-dependent (`p:%f`)
//! definitions. The six entries per family and body are, 1-based as in
//! the Fortran: `n:` nodal DOFs, `e:` edge DOFs, `f:` face DOFs, `d:` DG
//! DOFs, `b:` bubble DOFs and `p:` polynomial degree. Unset entries are
//! `-1` except `n:`, which defaults to `1`.

use crate::Error;

/// Element families 1..=10: point, line, tri, quad, tetra, pyramid, prism,
/// brick, and the face-bubble pseudo families `tri_face` (9) and
/// `quad_face` (10).
pub const FAMILIES: usize = 10;

/// Entry of a `Def_Dofs` row (1-based kind as in the Fortran).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Kind {
    Nodal = 1,
    Edge = 2,
    Face = 3,
    Dg = 4,
    Bubble = 5,
    P = 6,
}

/// One `[n, e, f, d, b, p]` row.
pub type Row = [i32; 6];

/// `inDOFs(10, 6)`: the mesh-wide maxima over all solvers and bodies that
/// `NonNodalElements` sizes the mesh from.
pub type MaxDefs = [Row; FAMILIES];

/// The unset `inDOFs` of `NonNodalElements`/`LoadModel`: `n:1`, rest `-1`.
pub fn initial_max_defs() -> MaxDefs {
    [[1, -1, -1, -1, -1, -1]; FAMILIES]
}

/// `Solver % Def_Dofs` of one solver.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ElementDefs {
    n_bodies: usize,
    /// `rows[(family - 1) * n_bodies + body - 1]`.
    rows: Vec<Row>,
}

/// Families a definition piece targets, by its leading family name
/// (`SEQL`, a prefix match): all eight basic families when none is given.
fn target_families(piece: &str) -> std::ops::RangeInclusive<usize> {
    const NAMES: [(&str, usize); 10] = [
        ("point", 1),
        ("line", 2),
        ("tri", 3),
        ("quad", 4),
        ("tetra", 5),
        ("pyramid", 6),
        ("prism", 7),
        ("brick", 8),
        ("tri_face", 9),
        ("quad_face", 10),
    ];
    // Later matches override earlier ones, as the Fortran's sequence of IFs
    // does ("tri_face" also starts with "tri").
    let mut ind = 1..=8;
    for (name, f) in NAMES {
        if piece.starts_with(name) {
            ind = f..=f;
        }
    }
    ind
}

/// Fortran `INDEX(s, key)` followed by a list-directed `READ(s(j+2:), *) l`.
fn read_int_after(piece: &str, key: &str) -> Result<Option<i32>, Error> {
    let Some(j) = piece.find(key) else {
        return Ok(None);
    };
    let rest = piece[j + key.len()..].trim_start();
    let token: String = rest
        .chars()
        .take_while(|c| !(c.is_whitespace() || *c == ',' || *c == '/'))
        .collect();
    token.parse::<i32>().map(Some).map_err(|_| {
        Error::BadDefinition(format!("cannot read an integer after '{key}' in '{piece}'"))
    })
}

/// Pieces of a definition separated by `-` (the target-family flags), with a
/// single leading `-` dropped, as the `DO WHILE` loops of `LoadModel` and
/// `GetMaxDefs` split them.
pub fn pieces(def: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let mut rest = def.trim();
    loop {
        if let Some(r) = rest.strip_prefix('-') {
            rest = r;
        }
        match rest.find('-') {
            Some(k) => {
                out.push(&rest[..k]);
                rest = &rest[k + 1..];
            }
            None => {
                out.push(rest);
                return out;
            }
        }
    }
}

/// `p:%name`: the MATC function name (the rest of the piece, trimmed), if
/// this piece defines the degree with one.
pub fn matc_degree_function(piece: &str) -> Option<&str> {
    let j = piece.find("p:")?;
    let rest = piece[j + 2..].strip_prefix('%')?;
    Some(rest.trim())
}

/// `MeshDeps` of `NonNodalElements`: the first `p:` of the whole definition
/// is followed by `%`.
pub fn is_mesh_dependent(def: &str) -> bool {
    def.find("p:")
        .is_some_and(|j| def[j + 2..].starts_with('%'))
}

impl ElementDefs {
    /// The unset table: `n:1`, all else `-1`.
    pub fn new(n_bodies: usize) -> Self {
        Self {
            n_bodies,
            rows: vec![[1, -1, -1, -1, -1, -1]; FAMILIES * n_bodies],
        }
    }

    pub fn n_bodies(&self) -> usize {
        self.n_bodies
    }

    /// `Def_Dofs(family, body, kind)`, 1-based family and body.
    ///
    /// # Panics
    /// If `family` or `body` is out of range.
    pub fn get(&self, family: usize, body: usize, kind: Kind) -> i32 {
        self.row(family, body)[kind as usize - 1]
    }

    /// The whole `[n, e, f, d, b, p]` row of a family and body.
    pub fn row(&self, family: usize, body: usize) -> &Row {
        assert!(
            (1..=FAMILIES).contains(&family),
            "family {family} out of range"
        );
        assert!(
            (1..=self.n_bodies).contains(&body),
            "body {body} out of range"
        );
        &self.rows[(family - 1) * self.n_bodies + body - 1]
    }

    fn row_mut(&mut self, family: usize, body: usize) -> &mut Row {
        &mut self.rows[(family - 1) * self.n_bodies + body - 1]
    }

    /// `ANY(Def_Dofs(family, :, 6) > 0)`: the condition of
    /// `isActivePElement` for elements of this family.
    pub fn p_active(&self, family: usize) -> bool {
        (1..=self.n_bodies).any(|b| self.get(family, b, Kind::P) > 0)
    }

    /// `ANY(Def_Dofs(:, :, 6) >= 1)` (`InitialPermutation`).
    pub fn any_p(&self) -> bool {
        self.rows.iter().any(|r| r[5] >= 1)
    }

    /// `LoadModel` on the solver-section `Element` keyword (`None` when it
    /// is absent and the solver is not DG, which means `"n:1"`): `GetDefs`
    /// on each piece for every body. `max_defs` is updated like the local
    /// `Def_Dofs(:,:)` of `LoadModel` (solvers on the global mesh).
    ///
    /// A `p:%f` degree sets `0` here; [`crate::PMesh`] evaluates it per
    /// element.
    pub fn from_solver_element(
        def: Option<&str>,
        n_bodies: usize,
        dg: bool,
        max_defs: &mut MaxDefs,
    ) -> Result<Self, Error> {
        let mut defs = Self::new(n_bodies);
        let def = match def {
            Some(d) => d,
            None if dg => {
                for f in 1..=FAMILIES {
                    for b in 1..=n_bodies {
                        defs.row_mut(f, b)[3] = 0;
                    }
                    max_defs[f - 1][3] = max_defs[f - 1][3].max(0);
                }
                return Ok(defs);
            }
            None => "n:1",
        };
        for piece in pieces(def) {
            defs.get_defs(piece, max_defs, dg)?;
        }
        Ok(defs)
    }

    /// `GetDefs` for one piece.
    fn get_defs(&mut self, piece: &str, max_defs: &mut MaxDefs, dg: bool) -> Result<(), Error> {
        let ind = target_families(piece);
        let set =
            |defs: &mut Self, kind: usize, l: i32, max_defs: &mut MaxDefs, update_max: bool| {
                for f in ind.clone() {
                    for b in 1..=defs.n_bodies {
                        defs.row_mut(f, b)[kind] = l;
                    }
                    if update_max {
                        max_defs[f - 1][kind] = max_defs[f - 1][kind].max(l);
                    }
                }
            };
        for (kind, key) in [(0, "n:"), (1, "e:"), (2, "f:")] {
            if let Some(l) = read_int_after(piece, key)? {
                set(self, kind, l, max_defs, true);
            }
        }
        match read_int_after(piece, "d:")? {
            Some(l) => set(self, 3, if l == 0 { -1 } else { l }, max_defs, true),
            None if dg => set(self, 3, 0, max_defs, true),
            None => {}
        }
        if let Some(l) = read_int_after(piece, "b:")? {
            set(self, 4, l, max_defs, true);
        }
        if matc_degree_function(piece).is_some() {
            set(self, 5, 0, max_defs, false);
        } else if let Some(l) = read_int_after(piece, "p:")? {
            set(self, 5, l, max_defs, true);
        }
        Ok(())
    }

    /// `GetMaxDefs` for one body: applies every piece of `def` to `body`
    /// and raises `max_defs`. `p_of(function)` evaluates a `p:%function`
    /// degree (the MATC call `function(cx)` with `cx` = body id and element
    /// midpoint, already rounded with `NINT`); it is only called for such
    /// pieces.
    ///
    /// Note the Fortran raises `max_defs` for families 1..=8 (or all ten for
    /// `n:`) whatever the target family of the piece; kept.
    pub fn get_max_defs(
        &mut self,
        def: &str,
        body: usize,
        dg: bool,
        max_defs: &mut MaxDefs,
        p_of: &mut dyn FnMut(&str) -> Result<i32, Error>,
    ) -> Result<(), Error> {
        for piece in pieces(def) {
            let ind = target_families(piece);
            let set = |defs: &mut Self, kind: usize, l: i32| {
                for f in ind.clone() {
                    defs.row_mut(f, body)[kind] = l;
                }
            };
            if let Some(l) = read_int_after(piece, "n:")? {
                set(self, 0, l);
                for row in max_defs.iter_mut() {
                    row[0] = row[0].max(l);
                }
            }
            for (kind, key) in [(1, "e:"), (2, "f:")] {
                if let Some(l) = read_int_after(piece, key)? {
                    set(self, kind, l);
                    for row in &mut max_defs[..8] {
                        row[kind] = row[kind].max(l);
                    }
                }
            }
            match read_int_after(piece, "d:")? {
                Some(l) => {
                    let l = if l == 0 { -1 } else { l };
                    set(self, 3, l);
                    for row in &mut max_defs[..8] {
                        row[3] = row[3].max(l);
                    }
                }
                None if dg => {
                    set(self, 3, 0);
                    for row in &mut max_defs[..8] {
                        row[3] = row[3].max(0);
                    }
                }
                None => {}
            }
            if let Some(l) = read_int_after(piece, "b:")? {
                set(self, 4, l);
                for row in &mut max_defs[..8] {
                    row[4] = row[4].max(l);
                }
            }
            let p = match matc_degree_function(piece) {
                Some(f) => Some(p_of(f)?),
                None => read_int_after(piece, "p:")?,
            };
            if let Some(l) = p {
                for row in &mut max_defs[..8] {
                    row[5] = row[5].max(l);
                }
                set(self, 5, l);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plain_degree_applies_to_basic_families() {
        let mut max = initial_max_defs();
        let d = ElementDefs::from_solver_element(Some("p:2"), 2, false, &mut max).unwrap();
        for f in 1..=8 {
            assert_eq!(d.get(f, 1, Kind::P), 2);
            assert_eq!(d.get(f, 2, Kind::P), 2);
            assert_eq!(d.get(f, 2, Kind::Nodal), 1);
            assert_eq!(max[f - 1][5], 2);
        }
        assert_eq!(d.get(9, 1, Kind::P), -1);
        assert!(d.p_active(4));
        assert!(!d.p_active(9));
    }

    #[test]
    fn family_flags_and_pieces() {
        assert_eq!(pieces("-quad p:3 -tri p:2"), vec!["quad p:3 ", "tri p:2"]);
        let mut max = initial_max_defs();
        let d = ElementDefs::from_solver_element(Some("-quad p:3 -tri p:2"), 1, false, &mut max)
            .unwrap();
        assert_eq!(d.get(4, 1, Kind::P), 3);
        assert_eq!(d.get(3, 1, Kind::P), 2);
        assert_eq!(d.get(5, 1, Kind::P), -1);
        // "tri_face" also starts with "tri"; the later test wins.
        assert_eq!(target_families("tri_face b:1"), 9..=9);
    }

    #[test]
    fn matc_degree_is_deferred() {
        let def = "p:%bodywisep";
        assert!(is_mesh_dependent(def));
        assert_eq!(matc_degree_function(def), Some("bodywisep"));
        let mut max = initial_max_defs();
        let mut d = ElementDefs::from_solver_element(Some(def), 2, false, &mut max).unwrap();
        assert_eq!(d.get(4, 1, Kind::P), 0);
        assert_eq!(max[3][5], -1);
        d.get_max_defs(def, 2, false, &mut max, &mut |f| {
            assert_eq!(f, "bodywisep");
            Ok(8)
        })
        .unwrap();
        assert_eq!(d.get(4, 2, Kind::P), 8);
        assert_eq!(d.get(4, 1, Kind::P), 0);
        assert_eq!(max[3][5], 8);
    }

    #[test]
    fn zero_dg_count_means_unset() {
        let mut max = initial_max_defs();
        let d = ElementDefs::from_solver_element(Some("n:1 b:3 d:0"), 1, false, &mut max).unwrap();
        assert_eq!(d.get(4, 1, Kind::Dg), -1);
        assert_eq!(d.get(4, 1, Kind::Bubble), 3);
        assert!(read_int_after("p:x", "p:").is_err());
    }
}
