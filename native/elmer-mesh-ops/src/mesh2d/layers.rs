//! The two generating layers of Mesh2D: `VoronoiVertex` and `VoronoiSegment`
//! (the `MovingFront` layer type), with `VSVertex`.

use super::geom::*;
use super::tri::*;

const SQRT3: f64 = 1.7320508075688772935274463415059;

/// `VoronoiVertex::generate`: insert the circumcentre of the worst triangle
/// until all are small enough.
pub fn generate_vertex(d: &mut Delaunay, nodes: &mut Nodes, bg: &Delaunay) {
    let vs = d.vs;
    let mut actives = Pq::new();
    let all = d.all_vertices.clone();
    for &v in &all {
        if !(d.verts[v].deleted || d.verts[v].is_external()) {
            let r = bg.interpolate(nodes, d.verts[v].vx, d.verts[v].vy);
            if !d.verts[v].right_size(r / SQRT3) {
                actives.insert(&mut d.verts, vs, v);
            }
        }
    }
    while actives.size() > 0 {
        let vtx = actives.first(&mut d.verts, vs).expect("active vertex");
        let nd = nodes.mesh_node(d.verts[vtx].vx, d.verts[vtx].vy);
        if d.add_site(nodes, vtx, nd, false, false, false) {
            if d.deleted.len() == 1 {
                // A terrible hack in the C++ too.
                d.verts[vtx].deleted = false;
                d.deleted.clear();
            } else {
                let del = d.deleted.clone();
                actives.remove_relevants(&mut d.verts, vs, &del);
                let nv = d.new_vertices.clone();
                for &v in &nv {
                    d.verts[v].body = 1;
                    let r = bg.interpolate(nodes, d.verts[v].vx, d.verts[v].vy);
                    if !d.verts[v].right_size(r / SQRT3) {
                        actives.insert(&mut d.verts, vs, v);
                    }
                }
            }
        }
        d.recycle();
    }
}

fn nb(d: &Delaunay, v: usize, i: usize) -> usize {
    d.verts[v].vertices[i].expect("neighbour of an interior vertex")
}

/// `VSVertex::testIfActive`.
fn test_if_active(d: &mut Delaunay, actives: &mut Pq, v: usize) -> bool {
    let n: Vec<usize> = (0..3).map(|i| nb(d, v, i)).collect();
    if n.iter()
        .any(|&u| d.verts[u].is_external() || d.verts[u].vs_type == VsType::Accepted)
    {
        d.verts[v].vs_type = VsType::Active;
        actives.insert(&mut d.verts, true, v);
        return true;
    }
    false
}

/// `VSVertex::radiate`.
fn radiate(d: &mut Delaunay, actives: &mut Pq, v: usize) {
    for i in 0..3 {
        let u = nb(d, v, i);
        if !d.verts[u].is_external() && d.verts[u].vs_type == VsType::Waiting {
            d.verts[u].vs_type = VsType::Active;
            actives.insert(&mut d.verts, true, u);
        }
    }
}

/// `VSVertex::borderTest`.
fn border_test(d: &mut Delaunay, v: usize) -> bool {
    let all = (0..3).all(|i| {
        let u = nb(d, v, i);
        d.verts[u].vs_type == VsType::Accepted || d.verts[u].is_external()
    });
    if all {
        d.verts[v].vs_type = VsType::Accepted;
    }
    all
}

/// `VSVertex::computeNewCoordinates`: the position of the next node on the
/// front, `None` when the vertex has to wait.
fn compute_new_coordinates(
    d: &mut Delaunay,
    nodes: &Nodes,
    bg: &Delaunay,
    v: usize,
) -> Option<(f64, f64)> {
    let vert = &d.verts[v];
    let n = vert.nodes;
    let mut len = 2.0 * vert.radius;
    let mut pick = 0;
    let mut touched = false;
    for i in 0..3 {
        let u = &d.verts[nb(d, v, i)];
        if u.is_external() || u.vs_type == VsType::Accepted {
            let (a, b) = (&nodes.v[n[i]], &nodes.v[n[(i + 1) % 3]]);
            let dist = distance(a.x, a.y, b.x, b.y);
            if dist < len {
                len = dist;
                pick = i;
                touched = true;
            }
        }
    }
    if !touched {
        d.verts[v].vs_type = VsType::Waiting;
        return None;
    }
    let (vx, vy, radius) = (vert.vx, vert.vy, vert.radius);
    let a = &nodes.v[n[pick]];
    let b = &nodes.v[n[(pick + 1) % 3]];
    let mx = (a.x + b.x) / 2.0;
    let my = (a.y + b.y) / 2.0;
    let h = bg.interpolate(nodes, mx, my);
    let mut ny = b.x - a.x;
    let mut nx = b.y - a.y;
    nx = -nx;
    nx /= len;
    ny /= len;
    let (mut c, mut dd) = (h, h);
    let (mut lx, mut ly);
    loop {
        let aa = (len * len + c * c - dd * dd) / (2.0 * len);
        if aa > c || (len - aa) > dd {
            return Some((vx, vy));
        }
        let bb = (c * c - aa * aa).sqrt();
        lx = a.x + aa * ny + bb * nx;
        ly = a.y - aa * nx + bb * ny;
        let vx1 = (lx + a.x) / 2.0;
        let vy1 = (ly + a.y) / 2.0;
        let vx2 = (lx + b.x) / 2.0;
        let vy2 = (ly + b.y) / 2.0;
        let h1 = bg.interpolate(nodes, vx1, vy1);
        let h2 = bg.interpolate(nodes, vx2, vy2);
        let e1 = (distance(lx, ly, a.x, a.y) - h1).abs() / h1;
        let e2 = (distance(lx, ly, b.x, b.y) - h2).abs() / h2;
        c = 0.75 * c + 0.25 * h1;
        dd = 0.75 * dd + 0.25 * h2;
        let error = e1 + e2;
        if error <= 1e-2 {
            break;
        }
    }
    if c > radius && dd > radius {
        lx = vx;
        ly = vy;
    }
    Some((lx, ly))
}

/// `VoronoiSegment::generate`: advance a front of accepted triangles.
pub fn generate_segment(d: &mut Delaunay, nodes: &mut Nodes, bg: &Delaunay, tag: i32) {
    let mut actives = Pq::new();
    let all = d.all_vertices.clone();
    for &v in &all {
        if d.verts[v].deleted {
            continue;
        }
        if d.verts[v].is_external() {
            for i in 0..3 {
                let Some(u) = d.verts[v].vertices[i] else {
                    continue;
                };
                if !d.verts[u].is_external() && d.verts[u].vs_type == VsType::Waiting {
                    d.verts[u].vs_type = VsType::Active;
                    actives.insert(&mut d.verts, true, u);
                }
            }
        } else {
            let r = bg.interpolate(nodes, d.verts[v].vx, d.verts[v].vy) / SQRT3;
            if d.verts[v].right_size(r) {
                if d.verts[v].vs_type == VsType::Active {
                    actives.remove(&mut d.verts, true, v);
                }
                d.verts[v].vs_type = VsType::Accepted;
                for i in 0..3 {
                    let u = nb(d, v, i);
                    if !d.verts[u].is_external() && d.verts[u].vs_type == VsType::Waiting {
                        d.verts[u].vs_type = VsType::Active;
                        actives.insert(&mut d.verts, true, u);
                    }
                }
            }
        }
    }

    while actives.size() > 0 {
        let vtx = actives.first(&mut d.verts, true).expect("active vertex");
        if border_test(d, vtx) {
            continue;
        }
        let Some((nx, ny)) = compute_new_coordinates(d, nodes, bg, vtx) else {
            continue;
        };
        let nd = nodes.mesh_node(nx, ny);
        if d.add_site(nodes, vtx, nd, false, true, false) {
            if d.deleted.len() == 1 {
                d.verts[vtx].deleted = false;
                d.verts[vtx].vs_type = VsType::Accepted;
                radiate(d, &mut actives, vtx);
                d.deleted.clear();
            } else {
                let del = d.deleted.clone();
                actives.remove_relevants(&mut d.verts, true, &del);
                let nv = d.new_vertices.clone();
                for &v in &nv {
                    d.verts[v].body = tag;
                }
                for &v in &nv {
                    let r = bg.interpolate(nodes, d.verts[v].vx, d.verts[v].vy) / SQRT3;
                    if d.verts[v].right_size(r) {
                        d.verts[v].vs_type = VsType::Accepted;
                    } else if !test_if_active(d, &mut actives, v) {
                        d.verts[v].vs_type = VsType::Waiting;
                    }
                }
                for &v in &nv {
                    match d.verts[v].vs_type {
                        VsType::Waiting => {
                            test_if_active(d, &mut actives, v);
                        }
                        VsType::Accepted => radiate(d, &mut actives, v),
                        VsType::Active => {}
                    }
                }
            }
        }
        d.recycle();
    }
}
