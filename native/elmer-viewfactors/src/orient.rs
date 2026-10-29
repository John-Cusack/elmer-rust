//! Which way a radiating boundary element faces: `DirectedNormalVector` of
//! `fem/src/ViewFactors.F90`, for line elements in the (x, y) plane.

use crate::Error;

/// A parent element of a boundary element, as `GetParentInfo` sees it.
#[derive(Copy, Clone, Debug)]
pub struct ParentInfo {
    /// Body id of the parent.
    pub body: i32,
    /// Coordinates of the last parent node not on the boundary element.
    pub node: [f64; 3],
}

/// The sign of `DirectedNormalVector` for the line element from `a` to `b`:
/// positive when the element's normal, the tangent turned a quarter
/// anticlockwise (`NormalVector`), points away from the body it radiates
/// from. See [`surface_direction`].
///
/// `target_body` is `Radiation Target Body` (or `Normal Target Body`) of the
/// boundary condition, `None` when the keyword is absent; a non-positive
/// value means "the other side of the meshed body". `radiative` gives
/// whether the `left`/`right` parent body has `Radiative Body = True`.
pub fn line_direction(
    a: [f64; 3],
    b: [f64; 3],
    left: Option<ParentInfo>,
    right: Option<ParentInfo>,
    target_body: Option<i32>,
    radiative: (bool, bool),
) -> Result<f64, Error> {
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let len = (dx * dx + dy * dy).sqrt();
    let nrm = [-dy / len, dx / len, 0.0];
    let centre = [
        (a[0] + b[0]) / 2.0,
        (a[1] + b[1]) / 2.0,
        (a[2] + b[2]) / 2.0,
    ];
    surface_direction(nrm, centre, left, right, target_body, radiative)
}

/// `DirectedNormalVector` for an element with normal `nrm` (`NormalVector`)
/// and centre `centre`: the sign of the normal's direction relative to the
/// body the element radiates into.
pub fn surface_direction(
    nrm: [f64; 3],
    centre: [f64; 3],
    left: Option<ParentInfo>,
    right: Option<ParentInfo>,
    target_body: Option<i32>,
    radiative: (bool, bool),
) -> Result<f64, Error> {
    let (mut lnode, mut rnode) = (left.map(|p| p.node), right.map(|p| p.node));
    let (lbody, rbody) = (left.map_or(0, |p| p.body), right.map_or(0, |p| p.body));

    let mut rad_body = match target_body {
        Some(b) => b,
        None => match radiative {
            (true, true) => {
                return Err(Error::Invalid(
                    "\"Radiative Body\" defined on both sides".into(),
                ))
            }
            (true, false) => rbody,
            (false, true) => lbody,
            (false, false) => 0,
        },
    };
    if rad_body < 0 {
        rad_body = 0;
    }
    if rad_body > 0 && rad_body != rbody && rad_body != lbody {
        return Err(Error::Invalid(
            "inconsistent direction information (Radiation Target Body)".into(),
        ));
    }

    if lnode.is_none() || (rad_body > 0 && rad_body == rbody) {
        std::mem::swap(&mut lnode, &mut rnode);
    }
    let mut normal_in = if rad_body <= 0 { 1.0 } else { -1.0 };
    // With a parent on both sides the radiating body's node is the more
    // reliable one under mesh deformation.
    if lnode.is_some() && rnode.is_some() {
        lnode = rnode;
        normal_in = -normal_in;
    }
    let interior = lnode.ok_or_else(|| Error::Invalid("boundary element has no parent".into()))?;
    let away = [
        centre[0] - interior[0],
        centre[1] - interior[1],
        centre[2] - interior[2],
    ];
    Ok(normal_in * (nrm[0] * away[0] + nrm[1] * away[1] + nrm[2] * away[2]))
}
